#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

namespace py = pybind11;

void fill_exclusion_mask(py::array_t<bool> py_mask,
                         const py::array_t<double> py_delta,
                         const py::array_t<double> py_oe,
                         const py::array_t<double> py_atom_coords,
                         const py::array_t<double> py_vdw_radii,
                         const double max_vdw_radius,
                         const double probe_radius)
{
    auto mask = py_mask.mutable_unchecked<3>();
    auto delta = py_delta.unchecked<1>();
    auto atom_coords = py_atom_coords.unchecked<2>();
    auto vdw_radii = py_vdw_radii.unchecked<1>();
    auto oe = py_oe.unchecked<1>();

    double centered_atom_x, centered_atom_y, centered_atom_z;
    double exclusion_distance, exclusion_distance2;
    double euclidean_z, z_distance, euclidean_y, y_distance, euclidean_x, x_distance, distance2;
    int32_t current_x, current_y, current_z;
    int32_t atom_idx, z_region, y_region, x_region;
    int32_t grid_x_coord, grid_y_coord, grid_z_coord;

    double half_delta_x = delta(0) / 2.0,
           half_delta_y = delta(1) / 2.0,
           half_delta_z = delta(2) / 2.0;

    double max_eclusion_dist = max_vdw_radius + probe_radius;

    // maximum distance (in cells) we want to check from a cell
    int32_t x_limit = (max_eclusion_dist + delta(0) * 1.5) / delta(0),
            y_limit = (max_eclusion_dist + delta(1) * 1.5) / delta(1),
            z_limit = (max_eclusion_dist + delta(2) * 1.5) / delta(2);

    for (atom_idx = 0; atom_idx < atom_coords.shape(0); atom_idx++)
    {
        // exclusion distance. van der waals radius + probe_radius (probe or probe + skin)
        exclusion_distance = vdw_radii(atom_idx) + probe_radius;
        exclusion_distance2 = exclusion_distance * exclusion_distance;

        // get grid cell coords from the atom's euclidean coords
        grid_x_coord = (atom_coords(atom_idx, 0) + half_delta_x - oe(0)) / delta(0);
        grid_y_coord = (atom_coords(atom_idx, 1) + half_delta_y - oe(1)) / delta(1);
        grid_z_coord = (atom_coords(atom_idx, 2) + half_delta_z - oe(2)) / delta(2);

        for (z_region = -z_limit; z_region <= z_limit; z_region++)
        {
            current_z = z_region + grid_z_coord;
            // check if the cell falls out of the grid
            if (current_z <= 0 or current_z > mask.shape(2))
                continue;

            // discrete grid coords back to euclidean space to check distances
            euclidean_z = current_z * delta(2) + oe(2);
            z_distance = euclidean_z - atom_coords(atom_idx, 2);
            if (z_distance > exclusion_distance)
                continue;

            for (y_region = -y_limit; y_region <= y_limit; y_region++)
            {
                current_y = y_region + grid_y_coord;
                // check if the cell falls out of the grid
                if (current_y <= 0 or current_y > mask.shape(1))
                    continue;

                euclidean_y = current_y * delta(1) + oe(1);
                y_distance = euclidean_y - atom_coords(atom_idx, 1);
                if (y_distance > exclusion_distance)
                    continue;

                for (x_region = -x_limit; x_region <= x_limit; x_region++)
                {
                    current_x = x_region + grid_x_coord;
                    // check if the cell falls out of the grid
                    if (current_x <= 0 or current_x > mask.shape(0))
                        continue;

                    if (mask(current_x - 1, current_y - 1, current_z - 1) == true)
                        continue;

                    euclidean_x = current_x * delta(0) + oe(0);
                    x_distance = euclidean_x - atom_coords(atom_idx, 0);
                    distance2 = x_distance * x_distance + y_distance * y_distance + z_distance * z_distance;
                    if (distance2 > exclusion_distance2)
                        continue;
                    mask(current_x - 1, current_y - 1, current_z - 1) = true;
                }
            }
        }
    }
}

PYBIND11_MODULE(skinlib, m)
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
