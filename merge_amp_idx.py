import xarray as xr
import numpy as np
import os

# User flag to enable splitting of specific boxes
flag_double = False  # Set to False if you need only 27 boxes

# Define box edges (x and y)
x_edges = [300, 421, 541, 661, 781, 901, 1021, 1141, 1261, 1307]
y_edges = [0, 128, 255, 380]

# Input and output paths
indir = "/work/cmcc/ag15419/basin_modes_num/"
outfile = "/work/cmcc/ag15419/basin_modes_num/basin_modes_amp_med.nc"

# Open the base dataset (box 1)
ds_full = xr.open_dataset(os.path.join(indir, "basin_modes_amp_1.nc")).copy(deep=True)

# Replace unwanted variable
if "sossheig" in ds_full:
    ds_full = ds_full.drop_vars("sossheig")

# Fill the entire domain with NaNs (except nav_lon/lat)
for var in ds_full.data_vars:
    if var not in ["nav_lon", "nav_lat"]:
        ds_full[var][:] = np.nan

# Boxes to split if flag_double is True (format: index: [(y0a, y1a), (y0b, y1b)])
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

# Loop over all 27 boxes
for row in range(3):        # 0 to 2 (south to north)
    for col in range(9):    # 0 to 8 (west to east)
        base_idx = row * 9 + col + 1
        x0, x1 = x_edges[col], x_edges[col + 1]
        y0, y1 = y_edges[row], y_edges[row + 1]

        if flag_double and base_idx in split_boxes:
            # Handle split boxes
            for i, (ys, ye) in enumerate(split_boxes[base_idx]):
                idx = f"{base_idx}{chr(97+i)}"  # 5a, 5b, etc.
                fname = f"basin_modes_amp_{idx}.nc"
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
                        print(f"Skipping variable {var}: unexpected shape {var_data.shape}")
        else:
            # Normal case
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
                else:
                    print(f"Skipping variable {var}: unexpected shape {var_data.shape}")

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

# Fix metadata: rename standard_name and units for phase variables
for var in ds_full.data_vars:
    if hasattr(ds_full[var], "standard_name") and ds_full[var].standard_name == "Phase":
        ds_full[var].attrs["standard_name"] = "Period"
        ds_full[var].attrs["units"] = "hours"

# Save output
ds_full.to_netcdf(outfile)
print(f"Final file saved as {outfile}")
