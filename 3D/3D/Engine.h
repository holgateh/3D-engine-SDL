#pragma once
#include<string>
#include<SDL.h>
#include"Algebra.h"
#include<iostream>
#include<chrono>
#include<ctime>
#include<fstream>
#include<algorithm>
#include<list>
struct Colour
{
	int r, g, b, a;
};

struct Camera
{
	Vec3 Pos;
	Vec3 Direction;

	Vec3 getLeftDirection()
	{
		return Direction.cross(Vec3(0, 1, 0)).normalise();
	}

	Vec3 getRightDirection()
	{
		return Direction.cross(Vec3(0, -1, 0)).normalise();
	}
};
class Engine
{
// Instance variables
private:
	SDL_Renderer* renderer;
	SDL_Window* window;
	int width, height;
	mat4x4 matProj;
	Camera vCamera;

	void update(float deltaTime); // do any updates to the scene movement/removal of objects.
	void render(SDL_Renderer* renderer); // display the updates to the screen.
	void drawTriangle(SDL_Renderer* renderer, triangle &tri, Colour col);
	void drawFilledTriangle(SDL_Renderer* renderer, triangle& tri, Colour col);

// Methods
public:
	Engine(int width, int height, std::string title);

	

};

