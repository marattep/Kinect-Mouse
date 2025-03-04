# Use the official Ubuntu image as the base image
FROM ubuntu:latest

# Install necessary packages
RUN apt-get update && apt-get install -y \
    cmake \
    build-essential \
    libgl1-mesa-dev \
    libglu1-mesa-dev \
    freeglut3-dev \
    libncurses5-dev libncursesw5-dev \
    libxtst-dev \
    libgl1-mesa-dev libglu1-mesa-dev libx11-dev \
    libfreenect-dev \
    libpthread-stubs0-dev \
    gdb 
    #only for debugging

#    libopencv-dev
    ## Install OpenCV ##########
# Install minimal prerequisites (Ubuntu 18.04 as reference)
#sudo apt update && sudo apt install -y cmake g++ wget unzip
# # Download and unpack sources
# wget -O opencv.zip https://github.com/opencv/opencv/archive/4.x.zip
# unzip opencv.zip
# # Create build directory
# mkdir -p build && cd build
# # Configure
# cmake  ../opencv-4.x
# # Build
# cmake --build .

#sudo apt-get install libopenni-dev


# Set the working directory
WORKDIR /app

# Command to run your application
#CMD ["cmake", ".", "&&", "make", "&&", "./kmouse"]

