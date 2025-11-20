import xarray as xr
import numpy as np
import os
from area_ini import *

# Define box edges (x and y)
x_edges = [300, 356, 412, 468, 524, 580, 636, 692, 748, 804,
           860, 916, 972, 1028, 1084, 1140, 1196, 1252, 1306]
y_edges = [0, 63, 126, 189, 252, 315, 379]

# Input and output paths
indir = work_dir
outfile = indir + "basin_modes_amp_med.nc"

# Step 1: Scan all files and collect all unique variable names (excluding nav_lon/nav_lat/sossheig)
print("Scanning files to collect all variable names...")
varnames_all = set()
shape_sample = None

for row in range(6):
    for col in range(18):
        base_idx = row * 18 + col + 1
        fname = f"basin_modes_amp_{base_idx}.nc"
        path = os.path.join(indir, fname)
        if not os.path.exists(path):
            continue
        ds_tmp = xr.open_dataset(path)
        for var in ds_tmp.data_vars:
            if var not in ["nav_lon", "nav_lat", "sossheig"]:
                varnames_all.add(var)
                if shape_sample is None:
                    shape_sample = ds_tmp[var].shape

varnames_all = sorted(varnames_all)
print(f"Collected {len(varnames_all)} data variables: {varnames_all}")

# Step 2: Open a base file to get coordinates and dimensions
basefile = os.path.join(indir, "basin_modes_amp_1.nc")
ds_base = xr.open_dataset(basefile)

# Create the empty dataset with proper coordinates
ds_full = xr.Dataset(coords={
    "y": ds_base.coords["y"],
    "x": ds_base.coords["x"],
    "nav_lon": ds_base["nav_lon"],
    "nav_lat": ds_base["nav_lat"]
})

# Add coordinates
for coord in ["nav_lon", "nav_lat"]:
    ds_full[coord] = ds_base[coord].copy(deep=True)

# Add NaN-filled variables
for var in varnames_all:
    template_var = None
    # Find a file where this var exists
    for row in range(6):
        for col in range(18):
            base_idx = row * 18 + col + 1
            path = os.path.join(indir, f"basin_modes_amp_{base_idx}.nc")
            if not os.path.exists(path):
                continue
            ds_tmp = xr.open_dataset(path)
            if var in ds_tmp:
                template_var = ds_tmp[var]
                break
        if template_var is not None:
            break

    if template_var is None:
        print(f"WARNING: Variable {var} not found in any input files. Skipping.")
        continue

    dims = template_var.dims
    shape = tuple(ds_base[dim].size if dim in ds_base.dims else template_var.sizes[dim] for dim in dims)
    data = np.full(shape, np.nan)
    ds_full[var] = xr.DataArray(data, dims=dims, coords={dim: ds_base[dim] for dim in dims if dim in ds_base.dims})
    ds_full[var].attrs = template_var.attrs

# Step 3: Fill data into `ds_full` as in original code
for row in range(6):
    for col in range(18):
        base_idx = row * 18 + col + 1
        x0, x1 = x_edges[col], x_edges[col + 1]
        y0, y1 = y_edges[row], y_edges[row + 1]

        fname = f"basin_modes_amp_{base_idx}.nc"
        path = os.path.join(indir, fname)
        if not os.path.exists(path):
            print(f"Missing file: {path}")
            continue

        ds_box = xr.open_dataset(path)
        for var in ds_box.data_vars:
            if var in ["nav_lon", "nav_lat", "sossheig"]:
                continue
            if var not in ds_full:
                continue
            var_data = ds_box[var].values
            if var_data.ndim == 2:
                ds_full[var].values[y0:y1, x0:x1] = var_data[y0:y1, x0:x1]
            elif var_data.ndim == 3:
                ds_full[var].values[:, y0:y1, x0:x1] = var_data[:, y0:y1, x0:x1]

# Set all values in box 73, 74, 91, 92(Biscay Bay) to NaN
exclude_boxes = [73, 74, 91, 92]

def box_bounds(idx, x_edges, y_edges, ncols=18):
    idx0 = idx - 1
    r = idx0 // ncols
    c = idx0 % ncols
    x0, x1 = x_edges[c], x_edges[c+1]
    y0, y1 = y_edges[r], y_edges[r+1]
    return x0, x1, y0, y1

for b in exclude_boxes:
    x0, x1, y0, y1 = box_bounds(b, x_edges, y_edges)
    for var in ds_full.data_vars:
        if var in ["nav_lon", "nav_lat"]:
            continue
        data = ds_full[var]
        if {"y", "x"}.issubset(data.dims):
            if "time_counter" in data.dims:
                data.values[:, y0:y1, x0:x1] = np.nan
            else:
                data.values[y0:y1, x0:x1] = np.nan

# Fix metadata
for var in ds_full.data_vars:
    if hasattr(ds_full[var], "standard_name") and ds_full[var].standard_name == "Phase":
        ds_full[var].attrs["standard_name"] = "Period"
        ds_full[var].attrs["units"] = "hours"

# Save output
ds_full.to_netcdf(outfile)
print(f"Final file saved as {outfile}")
