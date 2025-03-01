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

# Set the working directory
WORKDIR /app

# Command to run your application
#CMD ["cmake", ".", "&&", "make", "&&", "./kmouse"]

