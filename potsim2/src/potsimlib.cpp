#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <iostream>
#include <fstream>

namespace py = pybind11;

void fill_exclusion_mask(py::array_t<bool> &mask_array,
                         const py::array_t<double> &grid_spacing,
                         const py::array_t<double> &origin,
                         const py::array_t<double> &atom_positions,
                         const py::array_t<double> &atomic_radii,
                         const double &max_atomic_radius,
                         const double &probe_radius)
{
    // Get mutable/const references to numpy arrays
    auto mask = mask_array.mutable_unchecked<3>();
    const auto spacing = grid_spacing.unchecked<1>();
    const auto atoms = atom_positions.unchecked<2>();
    const auto radii = atomic_radii.unchecked<1>();
    const auto grid_origin = origin.unchecked<1>();

    // Precompute half grid spacings
    const double half_spacing_x = spacing(0) * 0.5;
    const double half_spacing_y = spacing(1) * 0.5;
    const double half_spacing_z = spacing(2) * 0.5;

    // Maximum distance to check for exclusion
    const double max_check_distance = max_atomic_radius + probe_radius;

    // Calculate grid cell limits for checking neighbors
    const int32_t max_cells_x = static_cast<int32_t>((max_check_distance + spacing(0) * 1.5) / spacing(0));
    const int32_t max_cells_y = static_cast<int32_t>((max_check_distance + spacing(1) * 1.5) / spacing(1));
    const int32_t max_cells_z = static_cast<int32_t>((max_check_distance + spacing(2) * 1.5) / spacing(2));

    // Process each atom
    for (int32_t atom_idx = 0; atom_idx < atoms.shape(0); ++atom_idx)
    {
        // Calculate exclusion radius for current atom
        const double exclusion_radius = radii(atom_idx) + probe_radius;
        const double exclusion_radius_sq = exclusion_radius * exclusion_radius;

        // Convert atom position to grid coordinates
        const int32_t grid_x = static_cast<int32_t>((atoms(atom_idx, 0) + half_spacing_x - grid_origin(0)) / spacing(0));
        const int32_t grid_y = static_cast<int32_t>((atoms(atom_idx, 1) + half_spacing_y - grid_origin(1)) / spacing(1));
        const int32_t grid_z = static_cast<int32_t>((atoms(atom_idx, 2) + half_spacing_z - grid_origin(2)) / spacing(2));

        // Iterate over neighboring grid cells
        for (int32_t dz = -max_cells_z; dz <= max_cells_z; ++dz)
        {
            const int32_t z = grid_z + dz;
            if (z <= 0 || z > mask.shape(2))
                continue;

            const double world_z = z * spacing(2) + grid_origin(2);
            const double delta_z = world_z - atoms(atom_idx, 2);
            if (std::abs(delta_z) > exclusion_radius)
                continue;

            for (int32_t dy = -max_cells_y; dy <= max_cells_y; ++dy)
            {
                const int32_t y = grid_y + dy;
                if (y <= 0 || y > mask.shape(1))
                    continue;

                const double world_y = y * spacing(1) + grid_origin(1);
                const double delta_y = world_y - atoms(atom_idx, 1);
                if (std::abs(delta_y) > exclusion_radius)
                    continue;

                for (int32_t dx = -max_cells_x; dx <= max_cells_x; ++dx)
                {
                    const int32_t x = grid_x + dx;
                    if (x <= 0 || x > mask.shape(0))
                        continue;

                    // Skip if cell is already marked
                    if (mask(x - 1, y - 1, z - 1))
                        continue;

                    const double world_x = x * spacing(0) + grid_origin(0);
                    const double delta_x = world_x - atoms(atom_idx, 0);

                    // Check if point is within exclusion sphere
                    const double dist_sq = delta_x * delta_x + delta_y * delta_y + delta_z * delta_z;
                    if (dist_sq <= exclusion_radius_sq)
                    {
                        mask(x - 1, y - 1, z - 1) = true;
                    }
                }
            }
        }
    }
}

PYBIND11_MODULE(potsimlib, m)
{
    m.doc() = R"pbdoc(
    pybind11 exclusion mask module
    -----------------------------

    Provides functionality for creating molecular exclusion masks.
    )pbdoc";

    m.def("fill_exclusion_mask", &fill_exclusion_mask,
          py::call_guard<py::gil_scoped_release>(),
          R"pbdoc(
          Fill exclusion mask for a set of atoms.
          
          Creates a spherical exclusion mask around atoms based on their 
          van der Waals radii and a probe radius.
          )pbdoc");
}
