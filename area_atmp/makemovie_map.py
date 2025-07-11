from cv2 import cv2
import os

workdir='/work/cmcc/ag15419/basin_modes/plots_area/'

image_folder = workdir 
video_name = workdir+'ssh_GB_BF.mp4'
images = [img for img in sorted(os.listdir(image_folder)) if img.endswith("Med_Gibraltar.png") and img.startswith("ssh_")] 
frame = cv2.imread(os.path.join(image_folder, images[0]))
height, width, layers = frame.shape

fourcc = cv2.VideoWriter_fourcc(*'mp4v')
video = cv2.VideoWriter(video_name, fourcc, 5, (width, height))

for image in images:
    print ('Infile: ',image)
    video.write(cv2.imread(os.path.join(image_folder, image)))

print ('Outfile: ',video_name)
cv2.destroyAllWindows()
video.release()
