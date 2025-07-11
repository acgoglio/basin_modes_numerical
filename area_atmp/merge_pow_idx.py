import xarray as xr
import numpy as np
import os
from area_ini import *

# User flag to enable splitting of specific boxes
flag_double = False  # Set to False if you need only 27 boxes (standard case)

# Define box edges (x and y)
x_edges = [300, 421, 541, 661, 781, 901, 1021, 1141, 1261, 1307]
y_edges = [0, 128, 255, 380]

# Input and output paths
indir = work_dir
outfile = indir+'basin_modes_pow_med.nc'

# Boxes to split if flag_double is True
split_boxes = {
    5: [(0, 64), (64, 128)],
    6: [(0, 64), (64, 128)],
    7: [(0, 64), (64, 128)],
    8: [(0, 64), (64, 128)],
    11: [(128, 191), (191, 255)],
    12: [(128, 191), (191, 255)],
    13: [(128, 191), (191, 255)],
    14: [(128, 191), (191, 255)],
    15: [(128, 191), (191, 255)],
    16: [(128, 191), (191, 255)],
}

# Step 1: Scan all files and collect all unique variable names (excluding nav_lon/nav_lat/sossheig)
print("Scanning files to collect all variable names...")
varnames_all = set()
shape_sample = None

for row in range(3):
    for col in range(9):
        base_idx = row * 9 + col + 1
        if flag_double and base_idx in split_boxes:
            for i in range(2):
                idx = f"{base_idx}{chr(97+i)}"
                fname = f"basin_modes_pow_{idx}.nc"
                path = os.path.join(indir, fname)
                if not os.path.exists(path):
                    continue
                ds_tmp = xr.open_dataset(path)
                for var in ds_tmp.data_vars:
                    if var not in ["nav_lon", "nav_lat", "sossheig"]:
                        varnames_all.add(var)
                        if shape_sample is None:
                            shape_sample = ds_tmp[var].shape
        else:
            fname = f"basin_modes_pow_{base_idx}.nc"
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
basefile = os.path.join(indir, "basin_modes_pow_1.nc")
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
    for row in range(3):
        for col in range(9):
            base_idx = row * 9 + col + 1
            idx_names = [f"{base_idx}a", f"{base_idx}b"] if flag_double and base_idx in split_boxes else [str(base_idx)]
            for idx in idx_names:
                path = os.path.join(indir, f"basin_modes_pow_{idx}.nc")
                if not os.path.exists(path):
                    continue
                ds_tmp = xr.open_dataset(path)
                if var in ds_tmp:
                    template_var = ds_tmp[var]
                    break
            if template_var is not None:
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
for row in range(3):
    for col in range(9):
        base_idx = row * 9 + col + 1
        x0, x1 = x_edges[col], x_edges[col + 1]
        y0, y1 = y_edges[row], y_edges[row + 1]

        if flag_double and base_idx in split_boxes:
            for i, (ys, ye) in enumerate(split_boxes[base_idx]):
                idx = f"{base_idx}{chr(97+i)}"
                fname = f"basin_modes_pow_{idx}.nc"
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
                        ds_full[var].values[ys:ye, x0:x1] = var_data[ys:ye, x0:x1]
                    elif var_data.ndim == 3:
                        ds_full[var].values[:, ys:ye, x0:x1] = var_data[:, ys:ye, x0:x1]
        else:
            fname = f"basin_modes_pow_{base_idx}.nc"
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

# Set all values in box 19 (Biscay Bay) to NaN
x0, x1 = x_edges[0], x_edges[1]
y0, y1 = y_edges[2], y_edges[3]
for var in ds_full.data_vars:
    if var in ["nav_lon", "nav_lat"]:
        continue
    data = ds_full[var]
    if {"y", "x"}.issubset(data.dims):
        if "time_counter" in data.dims:
            ds_full[var].values[:, y0:y1, x0:x1] = np.nan
        else:
            ds_full[var].values[y0:y1, x0:x1] = np.nan

# Fix metadata
for var in ds_full.data_vars:
    if hasattr(ds_full[var], "standard_name") and ds_full[var].standard_name == "Phase":
        ds_full[var].attrs["standard_name"] = "Period"
        ds_full[var].attrs["units"] = "hours"
    if hasattr(ds_full[var], "standard_name") and ds_full[var].standard_name == "Amplitude":
        ds_full[var].attrs["standard_name"] = "Energy"
        ds_full[var].attrs["units"] = "m2*s2"

# Save output
ds_full.to_netcdf(outfile)
print(f"Final file saved as {outfile}")
