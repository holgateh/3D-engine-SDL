#ifndef ALGEBRA_H
#define ALGEBRA_H
#include<vector>
struct Vec3
{
	float x = 0.0f;
	float y = 0.0f;
	float z = 0.0f;
	float w = 1.0f;

	Vec3() : x(0), y(0), z(0) {}
	Vec3(float _x, float _y, float _z) : x(_x), y(_y), z(_z), w(0) {}
	Vec3(const Vec3& v) : x(v.x), y(v.y), z(v.z), w(v.w) {}

	float mag()
	{
		return sqrt(x * x + y * y + z * z);
	}

	float dot(const Vec3& v)
	{
		return x * v.x + y * v.y + z * v.z;
	}

	Vec3 perp()
	{
		return Vec3(-y, x, 0);
	}

	Vec3 cross(const Vec3& v)
	{
		return Vec3(y * v.z - z * v.y, z * v.x - v.z * x, x * v.y - v.x * y);
	}

	Vec3 normalise()
	{
		float r = 1 / mag();
		return Vec3(x * r, y * r, z * r);
	}

	Vec3 operator + (const Vec3& rhs)
	{
		return Vec3(this->x + rhs.x, this->y + rhs.y, this->z + rhs.z);
	}

	Vec3 operator - (const Vec3& rhs)
	{
		return Vec3(this->x - rhs.x, this->y - rhs.y, this->z - rhs.z);
	}

	Vec3 operator * (const float& rhs)
	{
		return Vec3(rhs * this->x, rhs * this->y, rhs * this->z);
	}


};

struct triangle
{
	Vec3 p[3];
	wchar_t sym;
	float colour;
	Vec3 norm()
	{
		Vec3 n = p[1] - p[0];
		n = n.cross(p[2] - p[0]);
		return n;
	}

};

struct mesh
{
	std::vector<triangle> tris;
	Vec3 Pos = Vec3(0, 0, 0);
};

struct mat4x4
{
	float m[4][4] = { 0 };
};

struct mat3x3
{
	float m[3][3] = { 0 };
};


#endif