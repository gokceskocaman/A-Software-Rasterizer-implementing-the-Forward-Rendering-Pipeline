#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <cstring>
#include <string>
#include <vector>
#include <cmath>
#include <iostream>

#include "tinyxml2.h"
#include "Triangle.h"
#include "Helpers.h"
#include "Scene.h"

using namespace tinyxml2;
using namespace std;

/*
	Parses XML file
*/
Scene::Scene(const char *xmlPath)
{
	const char *str;
	XMLDocument xmlDoc;
	XMLElement *xmlElement;

	xmlDoc.LoadFile(xmlPath);

	XMLNode *rootNode = xmlDoc.FirstChild();

	// read background color
	xmlElement = rootNode->FirstChildElement("BackgroundColor");
	str = xmlElement->GetText();
	sscanf(str, "%lf %lf %lf", &backgroundColor.r, &backgroundColor.g, &backgroundColor.b);

	// read culling
	xmlElement = rootNode->FirstChildElement("Culling");
	if (xmlElement != NULL)
	{
		str = xmlElement->GetText();

		if (strcmp(str, "enabled") == 0)
		{
			this->cullingEnabled = true;
		}
		else
		{
			this->cullingEnabled = false;
		}
	}

	// read cameras
	xmlElement = rootNode->FirstChildElement("Cameras");
	XMLElement *camElement = xmlElement->FirstChildElement("Camera");
	XMLElement *camFieldElement;
	while (camElement != NULL)
	{
		Camera *camera = new Camera();

		camElement->QueryIntAttribute("id", &camera->cameraId);

		// read projection type
		str = camElement->Attribute("type");

		if (strcmp(str, "orthographic") == 0)
		{
			camera->projectionType = ORTOGRAPHIC_PROJECTION;
		}
		else
		{
			camera->projectionType = PERSPECTIVE_PROJECTION;
		}

		camFieldElement = camElement->FirstChildElement("Position");
		str = camFieldElement->GetText();
		sscanf(str, "%lf %lf %lf", &camera->position.x, &camera->position.y, &camera->position.z);

		camFieldElement = camElement->FirstChildElement("Gaze");
		str = camFieldElement->GetText();
		sscanf(str, "%lf %lf %lf", &camera->gaze.x, &camera->gaze.y, &camera->gaze.z);

		camFieldElement = camElement->FirstChildElement("Up");
		str = camFieldElement->GetText();
		sscanf(str, "%lf %lf %lf", &camera->v.x, &camera->v.y, &camera->v.z);

		camera->gaze = normalizeVec3(camera->gaze);
		camera->u = crossProductVec3(camera->gaze, camera->v);
		camera->u = normalizeVec3(camera->u);

		camera->w = inverseVec3(camera->gaze);
		camera->v = crossProductVec3(camera->u, camera->gaze);
		camera->v = normalizeVec3(camera->v);

		camFieldElement = camElement->FirstChildElement("ImagePlane");
		str = camFieldElement->GetText();
		sscanf(str, "%lf %lf %lf %lf %lf %lf %d %d",
			   &camera->left, &camera->right, &camera->bottom, &camera->top,
			   &camera->near, &camera->far, &camera->horRes, &camera->verRes);

		camFieldElement = camElement->FirstChildElement("OutputName");
		str = camFieldElement->GetText();
		camera->outputFilename = string(str);

		this->cameras.push_back(camera);

		camElement = camElement->NextSiblingElement("Camera");
	}

	// read vertices
	xmlElement = rootNode->FirstChildElement("Vertices");
	XMLElement *vertexElement = xmlElement->FirstChildElement("Vertex");
	int vertexId = 1;

	while (vertexElement != NULL)
	{
		Vec3 *vertex = new Vec3();
		Color *color = new Color();

		vertex->colorId = vertexId;

		str = vertexElement->Attribute("position");
		sscanf(str, "%lf %lf %lf", &vertex->x, &vertex->y, &vertex->z);

		str = vertexElement->Attribute("color");
		sscanf(str, "%lf %lf %lf", &color->r, &color->g, &color->b);

		this->vertices.push_back(vertex);
		this->colorsOfVertices.push_back(color);

		vertexElement = vertexElement->NextSiblingElement("Vertex");

		vertexId++;
	}

	// read translations
	xmlElement = rootNode->FirstChildElement("Translations");
	XMLElement *translationElement = xmlElement->FirstChildElement("Translation");
	while (translationElement != NULL)
	{
		Translation *translation = new Translation();

		translationElement->QueryIntAttribute("id", &translation->translationId);

		str = translationElement->Attribute("value");
		sscanf(str, "%lf %lf %lf", &translation->tx, &translation->ty, &translation->tz);

		this->translations.push_back(translation);

		translationElement = translationElement->NextSiblingElement("Translation");
	}

	// read scalings
	xmlElement = rootNode->FirstChildElement("Scalings");
	XMLElement *scalingElement = xmlElement->FirstChildElement("Scaling");
	while (scalingElement != NULL)
	{
		Scaling *scaling = new Scaling();

		scalingElement->QueryIntAttribute("id", &scaling->scalingId);
		str = scalingElement->Attribute("value");
		sscanf(str, "%lf %lf %lf", &scaling->sx, &scaling->sy, &scaling->sz);

		this->scalings.push_back(scaling);

		scalingElement = scalingElement->NextSiblingElement("Scaling");
	}

	// read rotations
	xmlElement = rootNode->FirstChildElement("Rotations");
	XMLElement *rotationElement = xmlElement->FirstChildElement("Rotation");
	while (rotationElement != NULL)
	{
		Rotation *rotation = new Rotation();

		rotationElement->QueryIntAttribute("id", &rotation->rotationId);
		str = rotationElement->Attribute("value");
		sscanf(str, "%lf %lf %lf %lf", &rotation->angle, &rotation->ux, &rotation->uy, &rotation->uz);

		this->rotations.push_back(rotation);

		rotationElement = rotationElement->NextSiblingElement("Rotation");
	}

	// read meshes
	xmlElement = rootNode->FirstChildElement("Meshes");

	XMLElement *meshElement = xmlElement->FirstChildElement("Mesh");
	while (meshElement != NULL)
	{
		Mesh *mesh = new Mesh();

		meshElement->QueryIntAttribute("id", &mesh->meshId);

		// read projection type
		str = meshElement->Attribute("type");

		if (strcmp(str, "wireframe") == 0)
		{
			mesh->type = WIREFRAME_MESH;
		}
		else
		{
			mesh->type = SOLID_MESH;
		}

		// read mesh transformations
		XMLElement *meshTransformationsElement = meshElement->FirstChildElement("Transformations");
		XMLElement *meshTransformationElement = meshTransformationsElement->FirstChildElement("Transformation");

		while (meshTransformationElement != NULL)
		{
			char transformationType;
			int transformationId;

			str = meshTransformationElement->GetText();
			sscanf(str, "%c %d", &transformationType, &transformationId);

			mesh->transformationTypes.push_back(transformationType);
			mesh->transformationIds.push_back(transformationId);

			meshTransformationElement = meshTransformationElement->NextSiblingElement("Transformation");
		}

		mesh->numberOfTransformations = mesh->transformationIds.size();

		// read mesh faces
		char *row;
		char *cloneStr;
		int v1, v2, v3;
		XMLElement *meshFacesElement = meshElement->FirstChildElement("Faces");
		str = meshFacesElement->GetText();
		cloneStr = strdup(str);

		row = strtok(cloneStr, "\n");
		while (row != NULL)
		{
			int result = sscanf(row, "%d %d %d", &v1, &v2, &v3);

			if (result != EOF)
			{
				mesh->triangles.push_back(Triangle(v1, v2, v3));
			}
			row = strtok(NULL, "\n");
		}
		mesh->numberOfTriangles = mesh->triangles.size();
		this->meshes.push_back(mesh);

		meshElement = meshElement->NextSiblingElement("Mesh");
	}
}

void Scene::assignColorToPixel(int i, int j, Color c)
{
	this->image[i][j].r = c.r;
	this->image[i][j].g = c.g;
	this->image[i][j].b = c.b;
}

/*
	Initializes image with background color
*/
void Scene::initializeImage(Camera *camera)
{
	if (this->image.empty())
	{
		for (int i = 0; i < camera->horRes; i++)
		{
			vector<Color> rowOfColors;
			vector<double> rowOfDepths;

			for (int j = 0; j < camera->verRes; j++)
			{
				rowOfColors.push_back(this->backgroundColor);
				rowOfDepths.push_back(1.01);
			}

			this->image.push_back(rowOfColors);
			this->depth.push_back(rowOfDepths);
		}
	}
	else
	{
		for (int i = 0; i < camera->horRes; i++)
		{
			for (int j = 0; j < camera->verRes; j++)
			{
				assignColorToPixel(i, j, this->backgroundColor);
				this->depth[i][j] = 1.01;
				this->depth[i][j] = 1.01;
				this->depth[i][j] = 1.01;
			}
		}
	}
}

/*
	If given value is less than 0, converts value to 0.
	If given value is more than 255, converts value to 255.
	Otherwise returns value itself.
*/
int Scene::makeBetweenZeroAnd255(double value)
{
	if (value >= 255.0)
		return 255;
	if (value <= 0.0)
		return 0;
	return (int)(value);
}

/*
	Writes contents of image (Color**) into a PPM file.
*/
void Scene::writeImageToPPMFile(Camera *camera)
{
	ofstream fout;

	fout.open(camera->outputFilename.c_str());

	fout << "P3" << endl;
	fout << "# " << camera->outputFilename << endl;
	fout << camera->horRes << " " << camera->verRes << endl;
	fout << "255" << endl;

	for (int j = camera->verRes - 1; j >= 0; j--)
	{
		for (int i = 0; i < camera->horRes; i++)
		{
			fout << makeBetweenZeroAnd255(this->image[i][j].r) << " "
				 << makeBetweenZeroAnd255(this->image[i][j].g) << " "
				 << makeBetweenZeroAnd255(this->image[i][j].b) << " ";
		}
		fout << endl;
	}
	fout.close();
}

/*
	Converts PPM image in given path to PNG file, by calling ImageMagick's 'convert' command.
*/
void Scene::convertPPMToPNG(string ppmFileName)
{
	string command;

	// TODO: Change implementation if necessary.
	//degistirmeyi sakinnnnn unutma
	//command = "convert " + ppmFileName + " " + ppmFileName + ".png";
	//system(command.c_str());
}

/*
	Transformations, clipping, culling, rasterization are done here.
*/



void Scene::applyTransformationsToVertices(vector<Vec4 *> *vertices4)
{

	//vector<Vec4 *> *vertices4Helper = & this->vertices4;
	for (int i = 0; i < vertices.size(); i++)
	{
		double x = vertices[i]->x, y = vertices[i]->y, z = vertices[i]->z;
		int colorId = vertices[i]->colorId;

		Vec4 *homogeneousVector = new Vec4(x, y, z, 1, colorId);
		vertices4->push_back(homogeneousVector);
		//Vec4 homogeneousVector(x, y, z, 1, colorId);
		//vertices4->push_back(&homogeneousVector);
	}

	for(int i = 0; i < meshes.size(); i++){
		vector<int> distinctVertices;
		//find vertices from "vertices" vector
		for(int j = 0; j < meshes[i]->triangles.size(); j++){
			if(find(distinctVertices.begin(), distinctVertices.end(), meshes[i]->triangles[j].vertexIds[0]) == distinctVertices.end()){
				distinctVertices.push_back(meshes[i]->triangles[j].vertexIds[0]);
			}
			if(find(distinctVertices.begin(), distinctVertices.end(), meshes[i]->triangles[j].vertexIds[1]) == distinctVertices.end()){
				distinctVertices.push_back(meshes[i]->triangles[j].vertexIds[1]);
			}
			if(find(distinctVertices.begin(), distinctVertices.end(), meshes[i]->triangles[j].vertexIds[2]) == distinctVertices.end()){
				distinctVertices.push_back(meshes[i]->triangles[j].vertexIds[2]);
			}
		}

		Matrix4 transformationMatrix;
		transformationMatrix = getIdentityMatrix();
		
		for(int j = 0; j < meshes[i] -> numberOfTransformations ; j++ ){
			char transformationType = meshes[i]->transformationTypes[j];
			if(transformationType == 't'){
				Translation translation = *translations[meshes[i]->transformationIds[j] - 1];
				
				double matrixHelper[4][4] = {{1, 0, 0, translation.tx},
								  			 {0, 1, 0, translation.ty},
								  			 {0, 0, 1, translation.tz},
								 			 {0, 0, 0, 1}};
				Matrix4 translationMatrix = Matrix4(matrixHelper);
				transformationMatrix = multiplyMatrixWithMatrix(translationMatrix, transformationMatrix);
			}
			else if(transformationType == 's'){
				Matrix4 scalingMatrix;
				Scaling scaling = *scalings[meshes[i]->transformationIds[j] - 1];
				double matrixHelper[4][4] = {{scaling.sx, 0, 0, 0},
								  			 {0, scaling.sy, 0, 0},
								  			 {0, 0, scaling.sz, 0},
								 			 {0, 0, 0, 1}};
				scalingMatrix = Matrix4(matrixHelper);
				transformationMatrix = multiplyMatrixWithMatrix(scalingMatrix, transformationMatrix);
			}
			else if(transformationType == 'r'){
				Rotation rotation = *rotations[meshes[i]->transformationIds[j] - 1];
				Vec3 u = Vec3(rotation.ux, rotation.uy, rotation.uz);
				Vec3 v;
				if(u.x != 0){
					v.x = u.y;
					v.y = -u.x;
					v.z = 0;
				}
				else{
					v.x = 0;
					v.y = u.z;
					v.z = -u.y;
				}
				Vec3 w = crossProductVec3(u, v);
				u = normalizeVec3(u);
				v = normalizeVec3(v);
				w = normalizeVec3(w);
				double baseToU[4][4] = {{u.x, u.y, u.z, 0},
										{v.x, v.y, v.z, 0},
										{w.x, w.y, w.z, 0},
										{0, 0, 0, 1}};
				double uToBase[4][4] = {{u.x, v.x, w.x, 0},
										{u.y, v.y, w.y, 0},
										{u.z, v.z, w.z, 0},
										{0, 0, 0, 1}};
				Matrix4 baseToUMatrix = Matrix4(baseToU);
				Matrix4 uToBaseMatrix = Matrix4(uToBase);
				double rotationMatrixHelper[4][4] = {{1, 0, 0, 0},
													 {0,  cos(M_PI * rotation.angle / 180), -sin(M_PI * rotation.angle / 180), 0},
													 {0, sin(M_PI * rotation.angle / 180),  cos(M_PI * rotation.angle / 180), 0},
													 {0, 0, 0, 1}};
				Matrix4 rotationMatrixHelperMatrix = Matrix4(rotationMatrixHelper);
				Matrix4 rotationMatrix = multiplyMatrixWithMatrix(uToBaseMatrix, rotationMatrixHelperMatrix);
				rotationMatrix = multiplyMatrixWithMatrix(rotationMatrix, baseToUMatrix);
				transformationMatrix = multiplyMatrixWithMatrix(rotationMatrix, transformationMatrix); 
			}
		}

		
		for(int j = 0; j < distinctVertices.size(); j++){
			int currentVertexId = (distinctVertices[j] - 1);
			Vec4 helperMatrix = * (vertices4->at(currentVertexId));
			//std::cout << helperMatrix.x << " " << helperMatrix.y << " " << helperMatrix.z << endl;
			Vec4 calculateTransformation = multiplyMatrixWithVec4(transformationMatrix, helperMatrix);
			Vec4 * calculateTransformationPointer = new Vec4(calculateTransformation.x, calculateTransformation.y, calculateTransformation.z, calculateTransformation.t, calculateTransformation.colorId);
			vertices4->at(currentVertexId) = calculateTransformationPointer;
			//std::cout << vertices4->at(currentVertexId)->x << " " << vertices4->at(currentVertexId)->y << " " << vertices4->at(currentVertexId)->z << endl;
		}
		
	}

}

void Scene::forwardRenderingPipeline(Camera *camera)
{
	// TODO: Implement this function

	// Apply modeling transformations to vertices
	std::vector<Vec4*> vertices4;
	applyTransformationsToVertices(&vertices4);

	for(int i = 0; i < vertices4.size(); i++){
		std::cout << vertices4.at(i)->x << " " << vertices4.at(i)->y << " " << vertices4.at(i)->z << endl;
	}






}
