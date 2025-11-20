from PIL import Image

# Location coordinates and name:
#where_point='312_531_atmp_Gulf_Lion'
#where_point='245_720_atmp_Tyrrhenian_Sea'
where_point='361_740_atmp_North_Adriatic'
where_point='95_711_atmp_Gulf_Gabes'

# Apro png files
img1 = Image.open("/work/cmcc/ag15419/basin_modes/basin_modes_num_atmp_nof/point_f/pow_"+where_point+".png").convert("RGBA")
img2 = Image.open("/work/cmcc/ag15419/basin_modes/basin_modes_num_atmp_nof/point/pow_"+where_point+".png").convert("RGBA")

# Stessa dimensione
img2 = img2.resize(img1.size)

# Sovrapposizione con trasparenza (blend)
blended = Image.blend(img1, img2, alpha=0.7)

blended.save("f-nof_"+where_point+".png")
