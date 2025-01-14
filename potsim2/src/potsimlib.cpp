#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include <iostream>
#include <fstream>

namespace py = pybind11;

void fill_exclusion_mask(py::array_t<bool> &py_mask,
                         const py::array_t<double> &py_delta,
                         const py::array_t<double> &py_oe,
                         const py::array_t<double> &py_atom_coords,
                         const py::array_t<double> &py_vdw_radii,
                         const double &max_vdw_radius,
                         const double &probe_radius)
{
    auto mask = py_mask.mutable_unchecked<3>();
    const auto delta = py_delta.unchecked<1>();
    const auto atom_coords = py_atom_coords.unchecked<2>();
    const auto vdw_radii = py_vdw_radii.unchecked<1>();
    const auto oe = py_oe.unchecked<1>();

    double exclusion_dist, exclusion_dist2;
    double z, z_dist, y, y_dist, x, x_dist, dist2;
    int32_t i, j, k;
    int32_t atom_idx, k_region, j_region, i_region;
    int32_t atom_i, atom_j, atom_k;

    const double half_delta_x = delta(0) / 2.0,
                 half_delta_y = delta(1) / 2.0,
                 half_delta_z = delta(2) / 2.0;

    const double max_exclusion_dist = max_vdw_radius + probe_radius;

    // maximum dist (in cells) we want to check from a cell
    const int32_t i_limit = (max_exclusion_dist + delta(0) * 1.5) / delta(0),
                  j_limit = (max_exclusion_dist + delta(1) * 1.5) / delta(1),
                  k_limit = (max_exclusion_dist + delta(2) * 1.5) / delta(2);

    for (atom_idx = 0; atom_idx < atom_coords.shape(0); atom_idx++)
    {
        // exclusion dist. van der waals radius + probe_radius (probe or probe + skin)
        exclusion_dist = vdw_radii(atom_idx) + probe_radius;
        exclusion_dist2 = exclusion_dist * exclusion_dist;

        // get grid cell coords from the atom's euclidean coords
        atom_i = (atom_coords(atom_idx, 0) + half_delta_x - oe(0)) / delta(0);
        atom_j = (atom_coords(atom_idx, 1) + half_delta_y - oe(1)) / delta(1);
        atom_k = (atom_coords(atom_idx, 2) + half_delta_z - oe(2)) / delta(2);

        for (k_region = -k_limit; k_region <= k_limit; k_region++)
        {
            k = k_region + atom_k;
            // check if the cell falls out of the grid
            if (k <= 0 or k > mask.shape(2))
                continue;

            // discrete grid coords back to euclidean space to check dists
            z = k * delta(2) + oe(2);
            z_dist = z - atom_coords(atom_idx, 2);
            if (z_dist > exclusion_dist)
                continue;

            for (j_region = -j_limit; j_region <= j_limit; j_region++)
            {
                j = j_region + atom_j;
                // check if the cell falls out of the grid
                if (j <= 0 or j > mask.shape(1))
                    continue;

                y = j * delta(1) + oe(1);
                y_dist = y - atom_coords(atom_idx, 1);
                if (y_dist > exclusion_dist)
                    continue;

                for (i_region = -i_limit; i_region <= i_limit; i_region++)
                {
                    i = i_region + atom_i;
                    // check if the cell falls out of the grid
                    if (i <= 0 or i > mask.shape(0))
                        continue;

                    if (mask(i - 1, j - 1, k - 1) == true)
                        continue;

                    x = i * delta(0) + oe(0);
                    x_dist = x - atom_coords(atom_idx, 0);
                    dist2 = x_dist * x_dist + y_dist * y_dist + z_dist * z_dist;
                    if (dist2 > exclusion_dist2)
                        continue;
                    mask(i - 1, j - 1, k - 1) = true;
                }
            }
        }
    }
}

PYBIND11_MODULE(potsimlib, m)
{
    m.doc() = R"pbdoc(
    pybind11 skin plugin
    -----------------------

    .. currentmodule:: skin

    .. autosummary::
        :toctree: _generate

        fill_exclusion_mask
    )pbdoc";

    m.def("fill_exclusion_mask", &fill_exclusion_mask, py::call_guard<py::gil_scoped_release>(), R"pbdoc(
        Fill exclusion mask
        
        Fill an spherical mask for given set of coordinates and radius.
    )pbdoc");
}
