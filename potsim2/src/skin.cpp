#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

namespace py = pybind11;

void fill_exclusion_mask(py::array_t<bool> pymask,
                         py::array_t<double> pydelta,
                         py::array_t<double> pyoe,
                         py::array_t<double> pypdb_center,
                         py::array_t<double> pyatom_coords,
                         py::array_t<double> pyvdw_radii,
                         double max_vdw,
                         double probe_radius)
{

    auto mask = pymask.mutable_unchecked<3>();
    auto delta = pydelta.unchecked<1>();
    auto atom_coords = pyatom_coords.unchecked<2>();
    auto vdw_radii = pyvdw_radii.unchecked<1>();
    auto oe = pyoe.unchecked<1>();
    auto pdb_center = pypdb_center.unchecked<1>();

    double cax, cay, caz, scax, scay, scaz;
    double rexcl, rexcl2;
    double ez, zdist, ey, ydist, ex, xdist, dist2;
    int x, y, z;
    int atom_idx, rz, ry, rx;
    int agcx, agcy, agcz;

    double half_deltax = delta(0) / 2.0,
           half_deltay = delta(1) / 2.0,
           half_deltaz = delta(2) / 2.0;

    double region_ex_max = max_vdw + probe_radius;

    // maximum cell distance we want to check from a cell
    int limx = (region_ex_max + delta(0) * 1.5) / delta(0),
        limy = (region_ex_max + delta(1) * 1.5) / delta(1),
        limz = (region_ex_max + delta(2) * 1.5) / delta(2);

    for (atom_idx = 0; atom_idx < atom_coords.shape(0); atom_idx++)
    {

        // center pdb structure in euclidean (0, 0, 0)
        cax = atom_coords(atom_idx, 0) - pdb_center(0);
        cay = atom_coords(atom_idx, 1) - pdb_center(1);
        caz = atom_coords(atom_idx, 2) - pdb_center(2);

        // shifting done in PIPSA
        scax = cax + half_deltax;
        scay = cay + half_deltay;
        scaz = caz + half_deltaz;

        // exlude distance. vdw + probe_radius(probe or probe + skin)
        rexcl = vdw_radii(atom_idx) + probe_radius;
        rexcl2 = rexcl * rexcl;

        // map euclidean atom coords into discrete grid cells
        agcx = (scax - oe(0)) / delta(0);
        agcy = (scay - oe(1)) / delta(1);
        agcz = (scaz - oe(2)) / delta(2);

        for (rz = -limz; rz <= limz; rz++)
        {
            z = rz + agcz;
            // checks if the cell falls out of the grid
            if (z <= 0 or z > mask.shape(2))
                continue;

            // discrete grid coords back to euclidean space to check distances
            // half delta shifting not corrected
            ez = z * delta(2) + oe(2);
            zdist = ez - caz;
            if (zdist > rexcl)
                continue;

            for (ry = -limy; ry <= limy; ry++)
            {
                y = ry + agcy;
                if (y <= 0 or y > mask.shape(1))
                    continue;

                ey = y * delta(1) + oe(1);
                ydist = ey - cay;
                if (ydist > rexcl)
                    continue;

                for (rx = -limx; rx <= limx; rx++)
                {
                    x = rx + agcx;
                    if (x <= 0 or x > mask.shape(0))
                        continue;

                    if (mask(x - 1, y - 1, z - 1) == true)
                        continue;

                    ex = x * delta(0) + oe(0);
                    xdist = ex - cax;
                    dist2 = xdist * xdist + ydist * ydist + zdist * zdist;
                    if (dist2 > rexcl2)
                        continue;

                    mask(x - 1, y - 1, z - 1) = true;
                }
            }
        }
    }
}

PYBIND11_MODULE(skinlib, m)
{
    m.doc() = R"pbdoc(
    Pybind11 skin plugin
    -----------------------

    .. currentmodule:: skin

    .. autosummary::
        :toctree: _generate

        fill_exclusion_mask
    )pbdoc";

    m.def("fill_exclusion_mask", fill_exclusion_mask, py::call_guard<py::gil_scoped_release>(), R"pbdoc(
        fill exclusion mask for a given radius
        
        .
    )pbdoc");
}
