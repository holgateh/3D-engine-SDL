# sdl-software-rasterizer

A software rasterizer built for fun using the SDL library. This is more of a learning project as I wanted to get my head around how a rasterizer and various other graphics pipeline steps work. SDL function calls to draw lines is all that is used to rasterize (very inefficient, I know). I have yet to implement a z-buffer for depth-testing, so the current implementation uses a rather hacky approximation: calculate the mean of distance from each of the triangle vertices to the camera; use this mean to sort the triangles into a drawing queue (closest traingles will be at the front of the queue). Another note: I should probably go back and tidy some of the code up as this was one of my first C++ projects, and some things are done in a rather amateurish fashion

![alt text](https://i.imgur.com/YRC6AvY.png)
