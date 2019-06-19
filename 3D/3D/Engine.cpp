#include "Engine.h"


std::vector<mesh> loadedModels;
std::vector<mesh> worldModels;
float fTheta = 0.0f;
bool normals = true;
float fYaw = 0.0f;
float fPitch = 0.0f;
float radius = 10.0f;


Vec3 lightSource = { 0, -1, 1, };


Vec3 Vector_IntersectPlane(Vec3& plane_p, Vec3& plane_n, Vec3& lineStart, Vec3& lineEnd)
{
	plane_n = plane_n.normalise();
	float plane_d = -plane_n.dot(plane_p);
	float ad = lineStart.dot(plane_n);
	float bd = lineEnd.dot(plane_n);
	float t = (-plane_d - ad) / (bd - ad);
	Vec3 lineStartToEnd = lineEnd - lineStart;
	Vec3 lineToIntersect = lineStartToEnd * t;
	return lineStart+ lineToIntersect;
}

int Triangle_ClipAgainstPlane(Vec3 plane_p, Vec3 plane_n, triangle& in_tri, triangle& out_tri1, triangle& out_tri2)
{
	// Make sure plane normal is indeed normal
	plane_n = plane_n.normalise();

	// Return signed shortest distance from point to plane, plane normal must be normalised
	auto dist = [&](Vec3& p)
	{
		Vec3 n = p.normalise();
		return (plane_n.x * p.x + plane_n.y * p.y + plane_n.z * p.z - plane_n.dot(plane_p));
	};

	// Create two temporary storage arrays to classify points either side of plane
	// If distance sign is positive, point lies on "inside" of plane
	Vec3* inside_points[3];  int nInsidePointCount = 0;
	Vec3* outside_points[3]; int nOutsidePointCount = 0;

	// Get signed distance of each point in triangle to plane
	float d0 = dist(in_tri.p[0]);
	float d1 = dist(in_tri.p[1]);
	float d2 = dist(in_tri.p[2]);

	if (d0 >= 0) { inside_points[nInsidePointCount++] = &in_tri.p[0]; }
	else { outside_points[nOutsidePointCount++] = &in_tri.p[0]; }
	if (d1 >= 0) { inside_points[nInsidePointCount++] = &in_tri.p[1]; }
	else { outside_points[nOutsidePointCount++] = &in_tri.p[1]; }
	if (d2 >= 0) { inside_points[nInsidePointCount++] = &in_tri.p[2]; }
	else { outside_points[nOutsidePointCount++] = &in_tri.p[2]; }

	if (nInsidePointCount == 0)
	{
		return 0; 
	}

	if (nInsidePointCount == 3)
	{
		out_tri1 = in_tri;
		return 1; 
	}

	if (nInsidePointCount == 1 && nOutsidePointCount == 2)
	{
		out_tri1.colour = in_tri.colour;
		out_tri1.sym = in_tri.sym;

		out_tri1.p[0] = *inside_points[0];

		out_tri1.p[1] = Vector_IntersectPlane(plane_p, plane_n, *inside_points[0], *outside_points[0]);
		out_tri1.p[2] = Vector_IntersectPlane(plane_p, plane_n, *inside_points[0], *outside_points[1]);

		return 1; 
	}

	if (nInsidePointCount == 2 && nOutsidePointCount == 1)
	{
		out_tri1.colour = in_tri.colour;
		out_tri1.sym = in_tri.sym;

		out_tri2.colour = in_tri.colour;
		out_tri2.sym = in_tri.sym;

		out_tri1.p[0] = *inside_points[0];
		out_tri1.p[1] = *inside_points[1];
		out_tri1.p[2] = Vector_IntersectPlane(plane_p, plane_n, *inside_points[0], *outside_points[0]);


		out_tri2.p[0] = *inside_points[1];
		out_tri2.p[1] = out_tri1.p[2];
		out_tri2.p[2] = Vector_IntersectPlane(plane_p, plane_n, *inside_points[1], *outside_points[0]);

		return 2; 
	}
}


float angle(Vec3 u, Vec3 v)
{
	return acosf((u.dot(v) / (u.mag() * v.mag())));
}

Vec3 multiplyMatrixVector4x4(Vec3 i, mat4x4&m)
{
	Vec3 v = Vec3();
	v.x = i.x * m.m[0][0] + i.y * m.m[1][0] + i.z * m.m[2][0] + i.w * m.m[3][0];
	v.y = i.x * m.m[0][1] + i.y * m.m[1][1] + i.z * m.m[2][1] + i.w * m.m[3][1];
	v.z = i.x * m.m[0][2] + i.y * m.m[1][2] + i.z * m.m[2][2] + i.w * m.m[3][2];
	v.w = i.x * m.m[0][3] + i.y * m.m[1][3] + i.z * m.m[2][3] + i.w * m.m[3][3];

	return v;
}


mat4x4 matrixMultiplymatrix(mat4x4 m1, mat4x4 m2)
{
	mat4x4 matrix;
	for (int c = 0; c < 4; c++)
		for (int r = 0; r < 4; r++)
			matrix.m[r][c] = m1.m[r][0] * m2.m[0][c] + m1.m[r][1] * m2.m[1][c] + m1.m[r][2] * m2.m[2][c] + m1.m[r][3] * m2.m[3][c];
	return matrix;
}


mat4x4 makeProjectionMatrix(float fAspectRatio, float fFovDegrees, float fFar, float fNear)
{
	float fFovRad = 1.0f / tanf(fFovDegrees * 0.5f / 180.0f * 3.14159f);
	mat4x4 matrix;
	matrix.m[0][0] = fAspectRatio * fFovRad;
	matrix.m[1][1] = fFovRad;
	matrix.m[2][2] = fFar / (fFar - fNear);
	matrix.m[3][2] = (-fFar * fNear) / (fFar - fNear);
	matrix.m[2][3] = 1.0f;
	matrix.m[3][3] = 0.0f;
	return matrix;
}

mat4x4 makeIdentity()
{
	mat4x4 matrix;
	matrix.m[0][0] = 1.0f;
	matrix.m[1][1] = 1.0f;
	matrix.m[2][2] = 1.0f;
	matrix.m[3][3] = 1.0f;
	return matrix;
}

mat4x4 makeRotationMatrixX(float theta)
{
	mat4x4 matrix;
	matrix.m[0][0] = 1.0f;
	matrix.m[1][1] = cosf(theta);
	matrix.m[1][2] = sinf(theta);
	matrix.m[2][1] = -sinf(theta);
	matrix.m[2][2] = cosf(theta);
	matrix.m[3][3] = 1.0f;
	return matrix;
}

mat4x4 makeRotationMatrixY(float theta)
{
	mat4x4 matrix;
	matrix.m[0][0] = cosf(theta);
	matrix.m[0][2] = sinf(theta);
	matrix.m[2][0] = -sinf(theta);
	matrix.m[1][1] = 1.0f;
	matrix.m[2][2] = cosf(theta);
	matrix.m[3][3] = 1.0f;
	return matrix;
}

mat4x4 makeRotationMatrixZ(float theta)
{
	mat4x4 matrix;
	matrix.m[0][0] = cosf(theta);
	matrix.m[0][1] = sinf(theta);
	matrix.m[1][0] = -sinf(theta);
	matrix.m[1][1] = cosf(theta);
	matrix.m[2][2] = 1.0f;
	matrix.m[3][3] = 1.0f;
	return matrix;
}

mat4x4 makeTranslationMatrix(float dx, float dy, float dz)
{
	mat4x4 matrix;
	matrix.m[0][0] = 1.0f;
	matrix.m[1][1] = 1.0f;
	matrix.m[2][2] = 1.0f;
	matrix.m[3][3] = 1.0f;
	matrix.m[3][0] = dx;
	matrix.m[3][1] = dy;
	matrix.m[3][2] = dz;
	return matrix;
}

mat4x4 makeQuickInverse(mat4x4 &m)
{
	mat4x4 matrix;
	matrix.m[0][0] = m.m[0][0]; matrix.m[0][1] = m.m[1][0]; matrix.m[0][2] = m.m[2][0]; matrix.m[0][3] = 0.0f;
	matrix.m[1][0] = m.m[0][1]; matrix.m[1][1] = m.m[1][1]; matrix.m[1][2] = m.m[2][1]; matrix.m[1][3] = 0.0f;
	matrix.m[2][0] = m.m[0][2]; matrix.m[2][1] = m.m[1][2]; matrix.m[2][2] = m.m[2][2]; matrix.m[2][3] = 0.0f;
	matrix.m[3][0] = -(m.m[3][0] * matrix.m[0][0] + m.m[3][1] * matrix.m[1][0] + m.m[3][2] * matrix.m[2][0]);
	matrix.m[3][1] = -(m.m[3][0] * matrix.m[0][1] + m.m[3][1] * matrix.m[1][1] + m.m[3][2] * matrix.m[2][1]);
	matrix.m[3][2] = -(m.m[3][0] * matrix.m[0][2] + m.m[3][1] * matrix.m[1][2] + m.m[3][2] * matrix.m[2][2]);
	matrix.m[3][3] = 1.0f;
	return matrix;
}

mat4x4 makePointAtMatrix(Vec3& pos, Vec3& target, Vec3& up)
{
	// Calculate new forward direction
	Vec3 newForward = (target - pos).normalise();

	// Calculate new Up direction
	Vec3 a = newForward * up.dot(newForward);
	Vec3 newUp = (up - a).normalise();

	// New Right direction is easy, its just cross product
	Vec3 newRight = newForward.cross(newUp);

	// Construct Dimensioning and Translation Matrix	
	mat4x4 matrix;
	matrix.m[0][0] = newRight.x;	matrix.m[0][1] = newRight.y;	matrix.m[0][2] = newRight.z;	matrix.m[0][3] = 0.0f;
	matrix.m[1][0] = newUp.x;		matrix.m[1][1] = newUp.y;		matrix.m[1][2] = newUp.z;		matrix.m[1][3] = 0.0f;
	matrix.m[2][0] = newForward.x;	matrix.m[2][1] = newForward.y;	matrix.m[2][2] = newForward.z;	matrix.m[2][3] = 0.0f;
	matrix.m[3][0] = pos.x;			matrix.m[3][1] = pos.y;			matrix.m[3][2] = pos.z;			matrix.m[3][3] = 1.0f;
	return matrix;
}

mesh loadModel(std::string filePath)
{
	std::vector<Vec3> vertices;
	std::vector<triangle> tris;
	mesh model;
	
	std::string line;
	std::ifstream file(filePath);
	if (file.is_open())
	{
		while (std::getline(file, line))
		{
			//std::cout << line << std::endl;
			// Add vertex to vertex list
			if (line[0] == 'v' && line.length() > 2)
			{
				Vec3 vertex;
				int i = 2;
				std::string comp[3];
				for (int j = 0; j < 3; j++)
				{
					while (i < line.length() && line[i] != ' ')
					{
						comp[j] += line[i];
						i++;
					}
					i++;
				}

				vertex.x = stof(comp[0]);
				vertex.y = stof(comp[1]);
				vertex.z = stof(comp[2]);
				vertices.push_back(vertex);
			}
			else if (line[0] == 'f' && line.length() > 2)
			{
				triangle tri;
				int i = 2;
				int indices[3];
				std::string comp[3];
				for (int j = 0; j < 3; j++)
				{
					while (i < line.length() && line[i] != ' ')
					{
						comp[j] += line[i];
						i++;
					}
					i++;
				}

				indices[0] = stoi(comp[0]);
				indices[1] = stoi(comp[1]);
				indices[2] = stoi(comp[2]);
				tri.p[0] = vertices[indices[0] - 1];
				tri.p[1] = vertices[indices[1] - 1];
				tri.p[2] = vertices[indices[2] - 1];
				tris.push_back(tri);
			}
		}
		model.tris = tris;
		file.close();
	}
	else
	{
		std::cout << "Unable to open file at " << filePath << std::endl;
	}

	return model;

}



void Engine::drawTriangle(SDL_Renderer* renderer, triangle &tri, Colour col)
{
	SDL_SetRenderDrawColor(renderer, col.r, col.g, col.b, SDL_ALPHA_OPAQUE);
	SDL_Point p[4];
	for (int i = 0; i < 3; i++)
	{
		p[i].x = tri.p[i].x;
		p[i].y = tri.p[i].y;
	}
	p[3] = p[0];

	SDL_RenderDrawLines(renderer, p, 4);
}

void Engine::drawFilledTriangle(SDL_Renderer* renderer, triangle &tri, Colour col)
{
	SDL_SetRenderDrawColor(renderer, col.r, col.g, col.b, SDL_ALPHA_OPAQUE);
	// Sort points of the triangle from min to max x;
	SDL_Point p[3] = {{tri.p[0].x, tri.p[0].y}, {tri.p[1].x, tri.p[1].y}, {tri.p[2].x, tri.p[2].y}};
	for (int i = 1; i < 3; i++)
	{
		int j;
		for (j = i - 1; j >= 0 && p[j].x > p[j+1].x; j--)
		{
			SDL_Point temp = p[j];
			p[j] = p[j+1];
			p[j+1] = temp;
		}
	}
	//std::cout << p[0].x << " " << p[1].x << " " << p[2].x << std::endl;

	int range = p[2].x - p[0].x;

	for (int i = 0; i < range; i++)
	{
		float x1, x2, y1, y2;
		x1 = 0;
		x2 = 0;
		y1 = 0;
		y2 = 0;
		if (i < range - (p[2].x - p[1].x) && range != 0.0f && p[1].x - p[0].x != 0.0f)
		{
			x1 = p[0].x + i;
			y1 = p[0].y + i * (p[2].y - p[0].y) / range;
			x2 = p[0].x + i;
			y2 = p[0].y + i * (p[1].y - p[0].y) / (p[1].x - p[0].x);
		}
		else if(range != 0.0f && p[2].x - p[1].x != 0.0f)
		{
			x1 = p[0].x + i;
			y1 = p[0].y + i * (p[2].y - p[0].y) / range;
			x2 = p[0].x + i;
			y2 = p[1].y + (i - (range - (p[2].x - p[1].x))) * (p[2].y - p[1].y) / (p[2].x - p[1].x);
		}
		SDL_RenderDrawLine(renderer, x1, y1, x2, y2);
	}
}

void Engine::update(float deltaTime)
{
	//fTheta += deltaTime * 0.00001f;
	mat4x4 rotY = makeRotationMatrixY(deltaTime * 0.1f);

	worldModels[1].Pos =  multiplyMatrixVector4x4(worldModels[1].Pos, rotY);

	for (auto& tri : worldModels[1].tris)
	{
		for (int i = 0; i < 3; i++)
		{
			tri.p[i] = multiplyMatrixVector4x4(tri.p[i], rotY);
		}
	}

	

}

void Engine::render(SDL_Renderer *renderer)
{
	// Clear the window before drawing
	SDL_SetRenderDrawColor(renderer, 0, 0, 0, SDL_ALPHA_OPAQUE);
	SDL_RenderClear(renderer);

	mat4x4 matRotZ, matRotX;
	matRotX = makeRotationMatrixX(0);
	matRotZ = makeRotationMatrixZ(0 * 0.5f);

	mat4x4 matTrans;
	matTrans = makeTranslationMatrix(10.0f , 2, 5.0f);

	mat4x4 matWorld;
	matWorld = makeIdentity();
	matWorld = matrixMultiplymatrix(matRotZ, matRotX);
	matWorld = matrixMultiplymatrix(matWorld, matTrans);

	
	Vec3 vUp = { 0, 1, 0 };
	Vec3 vTarget = { 0, 0, 1 };;
	
	mat4x4 cameraRotX = makeRotationMatrixX(fPitch);
	mat4x4 cameraRotY = makeRotationMatrixY(fYaw);

	vCamera.Direction = multiplyMatrixVector4x4(vTarget, cameraRotX).normalise();
	vCamera.Direction = multiplyMatrixVector4x4(vCamera.Direction, cameraRotY).normalise();
	vTarget = vCamera.Pos + vCamera.Direction;
	mat4x4 matCamera = makePointAtMatrix(vCamera.Pos, vTarget, vUp);


	mat4x4 matView = makeQuickInverse(matCamera);
	matView = matrixMultiplymatrix(matView, makeRotationMatrixZ(M_PI));


	std::vector<triangle> tris;

	for (auto& models : worldModels)
	{
		tris.insert(tris.end(), models.tris.begin(), models.tris.end());
	}


	std::vector<triangle> vecTrisToRaster;
	for (auto &tri : tris)
	{
		triangle triTransformed, triProjected, triViewed;
		for (int i = 0; i < 3; i++)
		{
			triTransformed.p[i] = multiplyMatrixVector4x4(tri.p[i], matWorld);
		}

		// Compute the normal direction to the plane containing the triangle.
		Vec3 norm = triTransformed.norm().normalise();
		// Only draw the triangle if the triangle is within our fov.
		if (norm.dot(vCamera.Direction) < 0.0f)
		{
			// Calculate illumation of the triangle.
			triViewed.colour = norm.dot(lightSource.normalise());
			if (triViewed.colour < 0)
				triViewed.colour = 0;

			// Scale into view;
			for (int i = 0; i < 3; i++)
			{
				triViewed.p[i] = multiplyMatrixVector4x4(triTransformed.p[i], matView);
			}
			int nClippedTriangles = 0;
			triangle clipped[2];
			nClippedTriangles = Triangle_ClipAgainstPlane({ 0.0f, 0.0f, 0.1f }, { 0.0f, 0.0f, 1.0f }, triViewed, clipped[0], clipped[1]);

			// We may end up with multiple triangles form the clip, so project as
			// required
			for (int n = 0; n < nClippedTriangles; n++)
			{
				for (int i = 0; i < 3; i++)
				{
					triProjected.p[i] = multiplyMatrixVector4x4(clipped[n].p[i],matProj);
					triProjected.p[i] = triProjected.p[i] * (1.0f / (triProjected.p[i].w));
					triProjected.p[i] = triProjected.p[i] + Vec3(1, 1, 0);
					triProjected.p[i].x *= 0.5f * Engine::width;
					triProjected.p[i].y *= 0.5f * Engine::height;
				}
				triProjected.colour = clipped[n].colour;
				triProjected.sym = clipped[n].sym;
				// Store triangle for sorting
				vecTrisToRaster.push_back(triProjected);
			}
		}

	}

	std::sort(vecTrisToRaster.begin(), vecTrisToRaster.end(), [](triangle& t1, triangle& t2)
		{
			float z1 = (t1.p[0].z + t1.p[1].z + t1.p[2].z) / 3.0f;
			float z2 = (t2.p[0].z + t2.p[1].z + t2.p[2].z) / 3.0f;
			return z1 > z2;
		});

		
	for (auto& triToRaster : vecTrisToRaster)
	{
		triangle clipped[2];
		std::list<triangle> listTriangles;

		// Add initial triangle
		listTriangles.push_back(triToRaster);
		int nNewTriangles = 1;

		for (int p = 0; p < 4; p++)
		{
			int nTrisToAdd = 0;
			while (nNewTriangles > 0)
			{
				// Take triangle from front of queue
				triangle test = listTriangles.front();
				listTriangles.pop_front();
				nNewTriangles--;

				// Clip it against a plane. We only need to test each 
				// subsequent plane, against subsequent new triangles
				// as all triangles after a plane clip are guaranteed
				// to lie on the inside of the plane. I like how this
				// comment is almost completely and utterly justified
				switch (p)
				{
				case 0:	nTrisToAdd = Triangle_ClipAgainstPlane({ 0.0f, 0.0f, 0.0f }, { 0.0f, 1.0f, 0.0f }, test, clipped[0], clipped[1]); break;
				case 1:	nTrisToAdd = Triangle_ClipAgainstPlane({ 0.0f, (float)Engine::height - 1, 0.0f }, { 0.0f, -1.0f, 0.0f }, test, clipped[0], clipped[1]); break;
				case 2:	nTrisToAdd = Triangle_ClipAgainstPlane({ 0.0f, 0.0f, 0.0f }, { 1.0f, 0.0f, 0.0f }, test, clipped[0], clipped[1]); break;
				case 3:	nTrisToAdd = Triangle_ClipAgainstPlane({ (float)Engine::width - 1, 0.0f, 0.0f }, { -1.0f, 0.0f, 0.0f }, test, clipped[0], clipped[1]); break;
				}

				// Clipping may yield a variable number of triangles, so
				// add these new ones to the back of the queue for subsequent
				// clipping against next planes
				for (int w = 0; w < nTrisToAdd; w++)
					listTriangles.push_back(clipped[w]);
			}
			nNewTriangles = listTriangles.size();
		}

		// Draw the final triangles
		for (auto& t : listTriangles)
		{
			drawFilledTriangle(renderer, t, { (int)(t.colour * 255), (int)(t.colour * 255), (int)(t.colour * 255), 0 });
			drawTriangle(renderer, t, { 0, 0, 0, 0 });
		}

	}
	


	//SDL_RenderDrawLine(renderer, 0, 0, Engine::width, Engine::height);
	SDL_RenderPresent(renderer);
}

Engine::Engine(int width, int height, std::string title)
{
	Engine::height = height;
	Engine::width = width;

	vCamera.Pos = Vec3(0, 0, 0);
	vCamera.Direction = Vec3(0, 0, 1);

	SDL_Point mousePos = { -1, -1 };

	//Load model

	std::vector<std::string> objFiles = { "Objects/terrian.obj", "Objects/earth.obj" , "Objects/monkey.obj"};
	for (auto &object : objFiles)
	{
		loadedModels.push_back(loadModel(object));
		worldModels.push_back(loadedModels[loadedModels.size() - 1]);

	}

	worldModels[1].Pos = Vec3(radius * sinf(fTheta), 0, radius * cosf(fTheta));

	mat4x4 transMatrix = makeTranslationMatrix(radius * sinf(fTheta), 0, radius * cosf(fTheta));

	for (auto& tri : worldModels[1].tris)
	{
		for (int i = 0; i < 3; i++)
		{
			tri.p[i] = multiplyMatrixVector4x4(tri.p[i], transMatrix);
		}
	}
	transMatrix = makeTranslationMatrix(0, -10.0f, 0);

	for (auto& tri : worldModels[0].tris)
	{
		for (int i = 0; i < 3; i++)
		{
			tri.p[i] = multiplyMatrixVector4x4(tri.p[i], transMatrix);
		}
	}


	// Initialise Projection matrix
	float fNear = 0.1f;
	float fFar = 1000.0f;
	float fFov = 90.0f;
	float fAspectRatio = height / (float)width;

	matProj = makeProjectionMatrix(fAspectRatio, fFov, fFar, fNear);

	

	// Initialise video and create window and renderer
	if (SDL_Init(SDL_INIT_EVERYTHING) == 0)
	{
		std::cout << "SDL_VIDEO initialised\n" << std::endl;
		SDL_SetRelativeMouseMode(SDL_TRUE);
		// Create pointers for renderer and window
		SDL_Window* window = NULL;
		SDL_Renderer* renderer = NULL;
		auto currentTime = std::chrono::system_clock::now();
		float elapsedTime = 0.0f;

		if (SDL_CreateWindowAndRenderer(width, height, 0, &window, &renderer) == 0)
		{

			bool isRunning = true;
			SDL_Event event;
			while (isRunning)
			{

				while (SDL_PollEvent(&event) != 0)
				{
					if (event.type == SDL_QUIT)
					{
						isRunning = false;
					}
					else if (event.type == SDL_KEYDOWN)
					{
						switch (event.key.keysym.sym)
						{
						case SDLK_a:
							vCamera.Pos = vCamera.Pos + vCamera.getLeftDirection() * elapsedTime;
							break;
						case SDLK_d:
							vCamera.Pos = vCamera.Pos + vCamera.getRightDirection() * elapsedTime;
							break;
						case SDLK_w:
							vCamera.Pos = vCamera.Pos + vCamera.Direction * 2.0f * elapsedTime;
							break;
						case SDLK_s:
							vCamera.Pos = vCamera.Pos - vCamera.Direction * 2.0f * elapsedTime;
							break;
						case SDLK_r:
							vCamera.Pos = vCamera.Pos + Vec3(0,1,0) * elapsedTime;
							break;
						case SDLK_f:
							vCamera.Pos = vCamera.Pos - Vec3(0, 1, 0) * elapsedTime;
							break;

						}
					}
					else if (event.type == SDL_MOUSEMOTION)
					{
							fYaw -= elapsedTime * event.motion.xrel * 1/20.0f;
							fPitch+= elapsedTime * event.motion.yrel * 1 /20.0f;		
							// Set pitch limiter.
							if (fPitch < -M_PI / 2.0f)
								fPitch = -M_PI / 2.0f;
							else if (fPitch > M_PI / 2.0f)
								fPitch = M_PI / 2.0f;
					}
				}

				//std::cout << "x = " << vCamera.Pos.x << " y = " << vCamera.Pos.y << " z = " << vCamera.Pos.z << std::endl;

				currentTime = std::chrono::system_clock::now();
				update(elapsedTime);
				render(renderer);
				auto newTime = std::chrono::system_clock::now();
				elapsedTime = (newTime - currentTime).count() * pow(10,-6);

				//std::cout << elapsedTime << std::endl;
			
			}
			SDL_DestroyWindow(window);
			SDL_Quit();

			
		}
	}

}