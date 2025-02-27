# Use the official Debian/Ubuntu image as the base image
#FROM ubuntu:latest

# Install necessary packages
RUN apt-get update && apt-get install -y \
    cmake \
    build-essential \
    # libgl1-mesa-dev \
    # libglu1-mesa-dev \
    # freeglut3-dev \
    # libx11-dev \
    # libxtst-dev \
    # libncurses5-dev \
    # libusb-1.0-0-dev \
    libfreenect-dev \
    libpthread-stubs0-dev


# Set the working directory
WORKDIR /app

# Copy your CMakeLists.txt and source code to the container
COPY . .

# Run CMake to configure the project
RUN cmake .

# Build the project
RUN make

# Command to run your application
CMD ["./kmouse"]
