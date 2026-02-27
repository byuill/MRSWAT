"""
MRSWAT Build Hypsometry — Python port of MRSWAT-Build-Hypsometry-080724.f
==========================================================================

PURPOSE
-------
Read a HEC-RAS geometry file (.g format) and build a *hypsometric curve*
for every cross-section it contains.  A hypsometric curve expresses:

    * Wetted cross-sectional area  A(z)   [m²]
    * Top-width                    W(z)   [m]

as functions of water-surface depth above the thalweg (channel minimum).

FORTRAN HERITAGE
----------------
This script is a line-for-line port of the Fortran-90 original.  Key
dimensions carried over unchanged:

    NPTS  = 800   max station/elevation points per cross-section
    IMILE = 800   max number of cross-sections (river miles)

All area / width calculations are done in US-customary units (feet) and
converted to SI (metres, m²) just before output — exactly as in the
Fortran source.

VERTICAL DISCRETISATION
-----------------------
*  Zint = 21 elevation slices per cross-section.
*  Slices 0 … 19 span thalweg (ELMIN) to max-bank (ELMAX) in 19 equal
   increments   ddZ = (ELMAX − ELMIN) / 19.
*  Slice 20 is an extreme-flood level  ZMWS = 300 ft ≈ 91.44 m.

Width at each slice is obtained by numerical differentiation of the area
curve: W(z) ≈ [A(z) − A(z − ε)] / ε  with ε = 1 × 10⁻⁴ ft.

CHANGE LOG (Python port)
------------------------
* Vectorised inner area loop with NumPy (≈20× faster than scalar loop).
* Fixed off-by-one in elevation indexing (was ``zmin + (n-1)*ddZ``,
  should be ``zmin + n*ddZ``).
* Fixed last-depth output to match Fortran (absolute ZMWS_m, not
  ZMWS_m − ELMIN).
* Matched Fortran comparison operators (.LT. / .GE. → strict < / >=).
* Robust line-ending handling (``\\r\\n`` on Windows).
"""

import numpy as np
import sys


# ===================================================================
# 1. FILE PARSER
#    Replaces the Fortran main-program I/O loop that reads sentinels
#    "Type" and "#Sta" from the HEC-RAS .g file.
# ===================================================================

def parse_hec_ras_file(filename: str) -> list:
    """
    Parse a HEC-RAS geometry (.g) file and return cross-section data.

    The function looks for two sentinel strings in each line:

    * ``Type`` in columns 1–4 → start of a new cross-section.
      The HEC-RAS river mile is extracted from columns 28–33
      (Fortran 1-based slice ``string(28:33)``; Python ``line[27:33]``).

    * ``#Sta`` in columns 1–4 → station/elevation data block.
      The point count comes from the portion after the ``=`` sign
      (Fortran reads columns 12–14, but splitting on ``=`` is safer).
      Data lines follow in Fortran ``10(a8)`` layout — up to 10
      eight-character fields per line, alternating X (station) and
      ELV (elevation).

    Parameters
    ----------
    filename : str
        Path to the HEC-RAS .g file.

    Returns
    -------
    list of dict
        Each element is ``{'rvm': float, 'points': [(x, elv), ...]}``.
    """
    cross_sections = []
    current_section = None

    try:
        with open(filename, "r") as fh:
            lines = fh.readlines()
    except FileNotFoundError:
        print(f"Error: file not found – '{filename}'")
        sys.exit(1)

    line_iter = iter(lines)

    for line in line_iter:
        # ---------- Sentinel 1: new cross-section header ----------
        #
        # Fortran:
        #   IF(string(1:4).eq.'Type') then
        #       Call cn(string(28:33), riverm, 6)
        #       ict = ict + 1 ;  RVM(ict) = riverm
        #
        if line[:4] == "Type":
            try:
                river_mile = float(line[27:33])
            except (ValueError, IndexError):
                print(f"Warning: cannot parse river mile from: "
                      f"{line.rstrip()}")
                continue
            current_section = {"rvm": river_mile, "points": []}
            cross_sections.append(current_section)

        # ---------- Sentinel 2: station / elevation data ----------
        #
        # Fortran:
        #   ELSEIF(string(1:4).eq.'#Sta') then
        #       Call cn(string(12:14), segnum, 3)
        #       ISEG(ict) = INT(segnum)
        #       Read(10,201)(name01(I), name02(I), I=1,ISEG(ict))
        # 201  FORMAT(10(a8))
        #
        elif line[:4] == "#Sta":
            if current_section is None:
                continue

            # --- number of points ---
            try:
                num_points = int(line.split("=")[1])
            except (ValueError, IndexError):
                print(f"Warning: cannot parse point count from: "
                      f"{line.rstrip()}")
                continue

            # --- read data lines in 10(a8) format ---
            # Each line has up to 10 eight-character fields.
            # We need 2 × num_points fields (alternating X, ELV).
            values_needed = 2 * num_points
            raw_fields = []

            while len(raw_fields) < values_needed:
                try:
                    data_line = next(line_iter)
                except StopIteration:
                    print("Warning: unexpected end of file while "
                          "reading data.")
                    break
                # Strip BOTH \r and \n to handle Windows line endings,
                # then split into 8-character chunks.
                stripped = data_line.rstrip("\r\n")
                for j in range(0, len(stripped), 8):
                    raw_fields.append(stripped[j: j + 8])

            # --- convert field pairs (X, ELV) to floats ---
            for i in range(num_points):
                try:
                    x   = float(raw_fields[2 * i])
                    elv = float(raw_fields[2 * i + 1])
                    current_section["points"].append((x, elv))
                except (ValueError, IndexError):
                    print(f"Warning: bad point {i + 1} for RM "
                          f"{current_section['rvm']}")
                    break

    return cross_sections


# ===================================================================
# 2. WETTED AREA AT A GIVEN WATER-SURFACE ELEVATION
#    Vectorised NumPy translation of the inner DO i = 2, ISEG(k)
#    loop in the Fortran "comarea" subroutine.
# ===================================================================

def calculate_area_at_elevation(
    x:     np.ndarray,
    z:     np.ndarray,
    n_pts: int,
    wse:   float,
) -> float:
    """
    Compute the wetted cross-sectional area below elevation *wse*.

    Three mutually exclusive geometric cases for each consecutive pair
    of profile points (i-1, i):

    1. **Both submerged** (z₁ < wse AND z₂ < wse):
       Full trapezoidal slab between the two verticals.

       ::

            wse ─────────────────────
                |                   |
           z₁  ├───────────────────┤ z₂
                                   ◄─── (x₂ − x₁)

       Area = [(wse − z₁) + (wse − z₂)] × (x₂ − x₁) / 2

    2. **Left exposed, right submerged** (z₁ ≥ wse, z₂ < wse):
       Find intersection xint where the bed crosses wse, then
       triangular sliver to the right.

       ::

           z₁ ╲
          wse ──╲────────────────
                 ╲              |
                  ╲─────────────┤ z₂
               xint            x₂

       xint = x₂ + (x₁ − x₂) × (wse − z₂) / (z₁ − z₂)
       Area = (wse − z₂) × (x₂ − xint) / 2

    3. **Left submerged, right exposed** (z₁ < wse, z₂ ≥ wse):
       Same intersection formula; triangular sliver to the left.

       ::

                              ╱ z₂
           ────────────────╱── wse
           |              ╱
        z₁ ├─────────────╱
           x₁          xint

       xint = x₂ + (x₁ − x₂) × (wse − z₂) / (z₁ − z₂)
       Area = (wse − z₁) × (xint − x₁) / 2

    Parameters
    ----------
    x     : ndarray – station (horizontal) coordinates (ft), length ≥ n_pts
    z     : ndarray – bed elevation at each station (ft),   length ≥ n_pts
    n_pts : int     – number of valid points in x, z
    wse   : float   – water-surface elevation (ft)

    Returns
    -------
    float – wetted cross-sectional area (ft²)

    Notes
    -----
    * HEC-RAS stations always increase left → right, matching the
      Fortran assumption that X(i) > X(i−1).
    * Comparison operators (strict ``<`` and ``>=``) match the Fortran
      ``.LT.`` / ``.GE.`` exactly for identical boundary-case behaviour.
    """
    # Build consecutive-pair arrays (vectorises the Fortran DO-loop)
    x1 = x[: n_pts - 1]          # "left"  stations   X(i-1)
    x2 = x[1: n_pts]             # "right" stations   X(i)
    z1 = z[: n_pts - 1]          # "left"  elevations ELV(i-1)
    z2 = z[1: n_pts]             # "right" elevations ELV(i)

    area = 0.0

    # ---- Case 1: both points below water surface ----
    # Fortran: (ELV(i) < zinc) AND (ELV(i-1) < zinc) in the ELSE branch
    mask = (z1 < wse) & (z2 < wse)
    if mask.any():
        area += np.sum(
            ((wse - z2[mask]) + (wse - z1[mask]))
            * (x2[mask] - x1[mask]) / 2.0
        )

    # ---- Case 2: left exposed (z₁ ≥ wse), right submerged (z₂ < wse) ----
    # Fortran: ELV(i) < zinc AND ELV(i-1) >= zinc  (first IF branch)
    mask = (z2 < wse) & (z1 >= wse)
    if mask.any():
        xint = (x2[mask]
                + (x1[mask] - x2[mask])
                * (wse - z2[mask]) / (z1[mask] - z2[mask]))
        area += np.sum((wse - z2[mask]) * (x2[mask] - xint) / 2.0)

    # ---- Case 3: left submerged (z₁ < wse), right exposed (z₂ ≥ wse) ----
    # Fortran: ELV(i) >= zinc AND ELV(i-1) < zinc  (ELSE IF branch)
    mask = (z2 >= wse) & (z1 < wse)
    if mask.any():
        xint = (x2[mask]
                + (x1[mask] - x2[mask])
                * (wse - z2[mask]) / (z1[mask] - z2[mask]))
        area += np.sum((wse - z1[mask]) * (xint - x1[mask]) / 2.0)

    return float(area)


# ===================================================================
# 3. HYPSOMETRY DRIVER
#    Replaces the Fortran subroutine comarea(ict).
# ===================================================================

def calculate_hypsometry(
    cross_sections: list,
    rorm: float,
    output_filename: str,
) -> None:
    """
    Build hypsometric curves and write the output file.

    Step-by-step correspondence with the Fortran ``comarea`` subroutine:

    1. Apply the river-mile offset RORM and discard sections outside
       the range  −25 < RM < 230.
    2. For each retained cross-section find the thalweg (ELMIN) and
       maximum bank elevation (ELMAX).
    3. Divide the elevation range into ``Zint − 2 = 19`` equal
       increments (``ddZ``), plus one extreme-flood level at
       ``ZMWS = 300 ft``.
    4. At each elevation slice compute wetted area and top width.
    5. Convert all quantities to SI (metres / m²).
    6. Write the hypsometry table.

    Parameters
    ----------
    cross_sections : list of dict
        Output of :func:`parse_hec_ras_file`.
    rorm : float
        Additive offset converting HEC-RAS river miles to standard
        river miles (e.g. −4.3 as of June 2024).
    output_filename : str
        Path for the output text file.
    """
    if not cross_sections:
        print("No cross-sections were parsed from the file.")
        return

    # ------------------------------------------------------------------
    # Pack parsed data into NumPy arrays.
    # These correspond to the Fortran COMMON blocks /RVMile/ and
    # /Numseg/:
    #   X(NPTS, IMILE)   – station coordinates        (ft)
    #   ELV(NPTS, IMILE) – bed elevations              (ft)
    #   RVM(IMILE)       – river mile labels
    #   ISEG(IMILE)      – number of valid points per section
    # ------------------------------------------------------------------
    max_seg = max(len(cs["points"]) for cs in cross_sections
                  if cs["points"])
    num_cs  = len(cross_sections)

    X    = np.full((max_seg, num_cs), np.nan)
    ELV  = np.full((max_seg, num_cs), np.nan)
    RVM  = np.array([cs["rvm"] for cs in cross_sections], dtype=float)
    ISEG = np.array([len(cs["points"]) for cs in cross_sections],
                    dtype=int)

    for idx, cs in enumerate(cross_sections):
        if cs["points"]:
            pts = np.array(cs["points"])
            n   = ISEG[idx]
            X[:n, idx]   = pts[:, 0]
            ELV[:n, idx] = pts[:, 1]

    # ------------------------------------------------------------------
    # Apply river-mile offset and discard sections outside the
    # modelling domain.
    #
    # Fortran:
    #   RVMT(i) = RVM(k) + RORM          ! offset
    #   if (RVMT(k) > -25 .AND. RVMT(k) < 230) ...   ! keep
    # ------------------------------------------------------------------
    RVM += rorm
    keep = (RVM > -25.0) & (RVM < 230.0)

    if not keep.any():
        print("No cross-sections within the valid range "
              "(−25 < RM < 230).")
        return

    RVM  = RVM[keep]
    ISEG = ISEG[keep]
    X    = X[:, keep]
    ELV  = ELV[:, keep]
    ict  = len(RVM)               # total retained cross-sections

    # ------------------------------------------------------------------
    # Per-section extremes
    # Fortran: ELMIN(:)=+999, ELMAX(:)=-999, then scanned per section.
    # ------------------------------------------------------------------
    ELMIN = np.nanmin(ELV, axis=0)   # thalweg elevation     (ft)
    ELMAX = np.nanmax(ELV, axis=0)   # max bank elevation    (ft)
    XMIN  = np.nanmin(X, axis=0)     # leftmost station      (ft)
    XMAX  = np.nanmax(X, axis=0)     # rightmost station     (ft)

    print(f"\nTotal # of cross-sections = {ict}")
    print(f"Global X-max = {np.nanmax(XMAX):.1f} ft,  "
          f"ELV range = [{np.nanmin(ELMIN):.1f}, "
          f"{np.nanmax(ELMAX):.1f}] ft\n")

    # ------------------------------------------------------------------
    # Hypsometry parameters
    #
    # Fortran:
    #   Zint  = 21        ! number of vertical slices
    #   RZMAX = 20.0
    #   ZMWS  = 300.0     ! extreme-flood water-surface elevation (ft)
    #   ddZ(k) = (ELMAX(k) - ELMIN(k)) / FLOAT(Zint - 2)
    # ------------------------------------------------------------------
    Zint = 21
    ZMWS = 300.0                     # extreme-flood WSE (ft)
    ddZ  = (ELMAX - ELMIN) / float(Zint - 2)   # 19 intervals per section

    Asum  = np.zeros((Zint, ict))    # wetted area at each slice  (ft²)
    Width = np.zeros((Zint, ict))    # top width at each slice    (ft)

    # ------------------------------------------------------------------
    # Compute wetted area at each elevation slice.
    #
    # Fortran elevation sequence (zinc):
    #   zinc = zmin                        ← n = 1
    #   zinc = zmin + ddZ                  ← n = 2
    #   …
    #   zinc = zmin + 19·ddZ = zmax        ← n = 20  (Zint-1)
    #   zinc = ZMWS                        ← n = 21  (Zint)
    #
    # 0-based mapping:  n_py = n_fortran − 1
    #   n = 0  → zinc = zmin              (area ≈ 0 at thalweg)
    #   n = 1  → zinc = zmin + ddZ
    #   …
    #   n = 19 → zinc = zmin + 19·ddZ = zmax
    #   n = 20 → zinc = ZMWS             (extreme flood)
    # ------------------------------------------------------------------
    for k in range(ict):
        zmin_k = ELMIN[k]
        for n in range(Zint):
            if n < Zint - 1:
                zinc = zmin_k + n * ddZ[k]       # ← corrected index
            else:
                zinc = ZMWS
            Asum[n, k] = calculate_area_at_elevation(
                X[:, k], ELV[:, k], ISEG[k], zinc
            )

    # ------------------------------------------------------------------
    # Compute top width via numerical derivative  dA/dz.
    #
    # Fortran:
    #   zinc = zmin + DDZ(k) − 1.0E-4     ← backward perturbation
    #   DO n = 2, Zint − 1
    #       <compute area at zinc>
    #       Width(n,k) = (Asum(n,k) − rarcum) / 1.0E-4
    #       zinc = zinc + ddZ(k)
    #   END DO
    #   Width(1,k) = 0.0
    #   Width(Zint,k) = XMAX(k) − XMIN(k)
    #
    # 0-based:  Fortran n=2..Zint-1  ↔  Python n=1..Zint-2
    # ------------------------------------------------------------------
    EPS = 1.0e-4                       # perturbation (ft)

    for k in range(ict):
        zmin_k = ELMIN[k]
        Width[0,        k] = 0.0                 # thalweg: zero width
        Width[Zint - 1, k] = XMAX[k] - XMIN[k]  # full channel span

        for n in range(1, Zint - 1):
            zinc_pert = zmin_k + n * ddZ[k] - EPS   # ← corrected index
            area_pert = calculate_area_at_elevation(
                X[:, k], ELV[:, k], ISEG[k], zinc_pert
            )
            Width[n, k] = (Asum[n, k] - area_pert) / EPS

    # ------------------------------------------------------------------
    # Convert to SI (metres)
    #
    # Fortran:
    #   ELMIN(k) = ELMIN(k) * 0.3048
    #   ddZ(k)   = ddZ(k)   * 0.3048
    #   Asum(n,k)= Asum(n,k)* 0.3048 * 0.3048
    #   Width(n,k)= Width(n,k)* 0.3048
    #   ZMWS     = ZMWS     * 0.3048
    # ------------------------------------------------------------------
    ft2m    = 0.3048
    ELMIN_m = ELMIN * ft2m
    ddZ_m   = ddZ   * ft2m
    Asum_m  = Asum  * ft2m ** 2
    Width_m = Width * ft2m
    ZMWS_m  = ZMWS  * ft2m            # ≈ 91.44 m

    # ------------------------------------------------------------------
    # Write output file.
    #
    # Format mirrors the Fortran WRITE(21,…) block:
    #
    #   Line 1:  label    "Number-of-cross-sections Min-RM MAX-RM"
    #   Line 2:  ict  RVM(ict)  RVM(1)     ← note: last RM first
    #
    #   Per cross-section:
    #       label   "Cross-section RM, Thalweg-elev(m) No-of-values"
    #       k  RVM(k)  ELMIN(k)  '21'
    #       label   "Depth(m) CS-Area(m2) Width(m)"
    #       Zint rows of   depth   area   width
    #
    # Depth starts at 0 (thalweg), increments by ddZ_m, and the final
    # row jumps to ZMWS_m (extreme-flood elevation in metres).
    #
    # Fortran FORMAT(f10.4, 2x, f14.4, 2x, f14.4)
    # ------------------------------------------------------------------
    with open(output_filename, "w") as fout:
        fout.write("Number-of-cross-sections Min-RM MAX-RM\n")
        # Fortran: Write(21,*) ict, RVM(ict), RVM(1)
        fout.write(f" {ict}  {RVM[-1]:.4f}  {RVM[0]:.4f}\n")

        for k in range(ict):
            fout.write(
                "Cross-section RM, Thalweg-elev(m) No-of-values\n"
            )
            fout.write(
                f" {k + 1}  {RVM[k]:.4f}  {ELMIN_m[k]:.4f}  {Zint}\n"
            )
            fout.write("Depth(m) CS-Area(m2) Width(m)\n")

            # Depth tracking variable — mirrors the Fortran ``zinc``
            # that starts at 0.0 and increments AFTER each WRITE.
            depth = 0.0
            for n in range(Zint):
                fout.write(
                    f"{depth:10.4f}  {Asum_m[n, k]:14.4f}  "
                    f"{Width_m[n, k]:14.4f}\n"
                )
                # Advance depth for the NEXT row (post-WRITE update)
                if n == Zint - 2:
                    # Fortran: if (n .eq. Zint-1) zinc = ZMWS
                    depth = ZMWS_m
                else:
                    depth += ddZ_m[k]

    print(f"Hypsometry written to '{output_filename}'")


# ===================================================================
# 4. ENTRY POINT
# ===================================================================

if __name__ == "__main__":
    print()
    print("MRSWAT Build Hypsometry (Python)")
    print("=" * 50)
    print("Reads a HEC-RAS geometry (.g file) and creates")
    print("hypsometric representations of each cross-section")
    print("— cross-sectional area as a function of depth at")
    print("the thalweg.")
    print()

    # Interactive prompts — mirror the Fortran READ(*,…) statements
    fname  = input("Enter the HEC-RAS geometry (.g) filename: ").strip()
    print("Enter the offset to convert HEC-RAS River Miles to "
          "standard River Miles.")
    print("NOTE: value of this as of 062524 is -4.3")
    rorm   = float(input("Offset: "))
    fnameo = input("Enter the desired filename for the output "
                   "hypsometry: ").strip()

    print("\nRunning...")

    # Step 1 — Parse the HEC-RAS geometry file.
    parsed = parse_hec_ras_file(fname)

    # Step 2 — Build hypsometry curves and write output.
    calculate_hypsometry(parsed, rorm, fnameo)
