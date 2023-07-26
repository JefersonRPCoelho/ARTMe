/**
*This software is provided by the Tecgraf/PUC-Rio Institute "as is" and any express or implied
*warranties, including, but not limited to, the implied warranties of merchantability and fitness
*for a particular purpose are disclaimed. In no event shall the scalable software infrastructure
*project be liable for any direct, indirect, incidental, special, exemplary, or consequential
*damages (including, but not limited to, procurement of substitute goods or services; loss of
*use, data, or profits; or business interruption) however caused and on any theory of liability,
*whether in contract, strict liability, or tort (including negligence or otherwise) arising in
*any way out of the use of this software, even if advised of the possibility of such damage.
*/


/* 
 * File:   main.cpp
 * Author: jcoelho
 *
 * A simple test that refines all triangles inside a hexagon and over a sphere using the different
 * triangulation patterns in order to illustrate the syntax of the refinement lib.
 *
 */
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <set>
#include <string>
#include <fstream>
#include "CornerTable.h"
#include "CornerTableAdaptiveRefinement.h"

#ifndef M_PI
#define M_PI  3.14159265358979323846
#endif

using namespace std;



/**
 * Create a mesh of an regular hexagon.
 * @return - corner table with the mesh.
 */
CornerTable* createHexagonMesh( )
{
    //Allocate the geometry vector.
    double *coordinates = new double[3 * 7];

    //Allocate the triangulation vector.
    CornerType *triangleMesh = new CornerType[3 * 6];

    //Central vertex.
    double centerX = 0.0, centerY = 0.0;

    //Radius of the hexagon.
    double radius = 1.0;

    //Set the coordinates of the central vertex.
    coordinates[0] = centerX;
    coordinates[1] = centerY;
    coordinates[2] = 0.0;

    //Set the first point in x-axe direction.
    coordinates[3] = centerX + radius;
    coordinates[4] = centerY;
    coordinates[5] = 0.0;

    for (unsigned int i = 1; i < 6; i++)
    {
        //Compute the angle.
        double angle = i * M_PI / 3.0;

        //Set the new vertex coordinates.
        coordinates[3 * ( i + 1 ) + 0] = radius * cos( angle );
        coordinates[3 * ( i + 1 ) + 1] = radius * sin( angle );
        coordinates[3 * ( i + 1 ) + 2] = 0.0;

        //Set new triangle.
        triangleMesh[3 * ( i - 1 ) + 0] = 0;
        triangleMesh[3 * ( i - 1 ) + 1] = i;
        triangleMesh[3 * ( i - 1 ) + 2] = i + 1;
    }
    //Set last triangle.
    triangleMesh[3 * 5 + 0] = 0;
    triangleMesh[3 * 5 + 1] = 6;
    triangleMesh[3 * 5 + 2] = 1;

    //Allocate the corner table.
    CornerTable *surface = new CornerTable( triangleMesh, coordinates, 6, 7, 3 );

    //Free memory.
    delete []triangleMesh;
    delete []coordinates;

    return surface;
}



std::set<CornerType> getTrianglesToBeRefined( CornerTableAdaptiveRefinement* refinement )
{
    //Get the surface representation.
    CornerTable* surface = refinement->getSurface( );

    //A set of selected triangles to be refined.
    std::set<CornerType> selectedTriangles;

    //Get the number of triangles on current mesh.
    CornerType numberTriangles = surface->getNumTriangles( );

    //Get the triangle mesh indexes.
    const CornerType* triangleMesh = surface->getTriangleList( );

    //Get the triangle mesh coordinates.
    const double* coords = surface->getAttributes( );

    //Get the number of coordinates by points.
    int d = surface->getNumberAttributesByVertex( );

    //Select triangles to be refined.
    for (CornerType t = 0; t < numberTriangles; t++)
    {
        for (int i = 0; i < 3; i++)
        {
            //Get the point index.
            CornerType p = triangleMesh[3 * t + i];

            //Add triangle that has points in the interval y < -0.7 || y > 0.7
            if (fabs( coords[d * p + 1] ) > 0.7)
            {
                selectedTriangles.insert( t );
                break;
            }
        }
    }

    return selectedTriangles;
}



/**
 * Callback to compute the coordinates of the new vertices. In this case, the
 * callback set new coordinates in the unit sphere surface.
 * @param surface - current surface represented in the corner table data
 * structure.
 * @param corner - corner opposite to the edge where the new vertex will be created.
 * @param attributes - vector to be set with the new coordinates.
 */
static void computeNewVertexOverSphere( const CornerTable* surface, const CornerType corner,
                                        double* attributes )
{
    //Get the index of the vertices on edge.
    const CornerType p1 = surface->cornerToVertexIndex( surface->cornerNext( corner ) );
    const CornerType p2 = surface->cornerToVertexIndex( surface->cornerPrevious( corner ) );

    //Get the geometry vector.
    const double* G = surface->getAttributes( );

    //Get the number of coordinates by vector.
    unsigned int d = surface->getNumberAttributesByVertex( );

    //Compute the middle point and set on vector.
    for (unsigned int i = 0; i < d; i++)
    {
        attributes[i] = ( G[d * p1 + i] + G[d * p2 + i] ) / 2.0;
    }

    //As the sphere has its center in (0, 0, 0), the solution is just compute the unit vector to
    //the middle point.
    double norm = 0.0;
    for (unsigned int i = 0; i < d; i++)
    {
        norm += attributes[i] * attributes[i];
    }

    norm = sqrt( norm );

    for (unsigned int i = 0; i < d; i++)
    {
        attributes[i] /= norm;
    }
}



/**
 * Refine the top and the bottom of the hexagon/sphere with different refinement levels
 * as constraints.
 * @param refinement - object to perform the refinement.
 * @param r - distance from the center to include triangles.
 */
void refinementStep( CornerTableAdaptiveRefinement* refinement, CornerTableAdaptiveRefinement::Triangulation pattern, bool sphere )
{

    //A set of selected triangles to be refined.
    std::set<CornerType> selectedTriangles = getTrianglesToBeRefined( refinement );

    //Refine the triangle mesh with undo mechanism, the ARTMe triangulation and
    //our callback. The radius is set as 0 because it is used just in the
    //incremental refinement.
    if (sphere)
        refinement->meshAdaptiveRefinement( selectedTriangles, true, pattern, 3, computeNewVertexOverSphere );
    else
        refinement->meshAdaptiveRefinement( selectedTriangles, true, pattern, 3 );
}



CornerTable* readOFFFIle( const char* path )
{
    std::ifstream in( path );
    if (in.fail( ))
    {
        printf( "Error opening the file %s\n", path );
        return 0;
    }

    CornerType numPoints, numberTriangles;
    unsigned int numAttrib = 3;

    std::string trash;
    in >> trash;
    in >> numPoints >> numberTriangles >> trash;

    double* coords = new double[3 * numPoints];
    CornerType* triangles = new CornerType[3 * numberTriangles];

    //Write the coordinates.
    for (CornerType i = 0; i < numPoints; i++)
    {
        for (unsigned int c = 0; c < numAttrib; c++)
        {
            in >> coords[numAttrib * i + c];
        }
    }

    for (CornerType i = 0; i < numberTriangles; i++)
    {
        in >> trash >> triangles[3 * i + 0] >> triangles[3 * i + 1] >> triangles[3 * i + 2];
    }

    //Allocate the corner table.
    CornerTable *surface = new CornerTable( triangles, coords, numberTriangles, numPoints, 3 );

    //Free memory.
    delete []triangles;
    delete []coords;

    return surface;
}



void refinementExample( CornerTableAdaptiveRefinement::Triangulation pattern, std::string name, bool sphere )
{
    CornerTable* surface = 0;
    if (sphere)
        surface = readOFFFIle( "../input/icosaedro.off" );
    else
        surface = createHexagonMesh( );

    //Create the refinement object.
    CornerTableAdaptiveRefinement* refinement = new CornerTableAdaptiveRefinement( surface );

    //Write the original mesh.
    refinement->writeMeshOFFFIle( ( name + "original.off" ).c_str( ) );

    LevelType maxLevel = 5;
    for (LevelType l = 0; l < maxLevel; l++)
    {
        refinementStep( refinement, pattern, sphere );
    }

    //Write the resulting mesh.
    refinement->writeMeshOFFFIle( ( name + "refined.off" ).c_str( ) );

    //Perform all undo steps.
    while (!refinement->undoIsEmpty( ))
    {
        refinement->makeUndo( );
    }

    //Write the resulting mesh after all undo steps.
    refinement->writeMeshOFFFIle( ( name + "originalCopy.off" ).c_str( ) );

    //Perform all redo steps.
    while (!refinement->redoIsEmpty( ))
    {
        refinement->makeRedo( );
    }

    //Write the resulting mesh after all redo steps.
    refinement->writeMeshOFFFIle( ( name + "refinedCopy.off" ).c_str( ) );

    //Free memory.
    delete refinement;
    delete surface;
}



/*
 * 
 */
int main( int argc, char** argv )
{
    //Hexagon example with three different patterns.
    {
        refinementExample( CornerTableAdaptiveRefinement::ARTMe, "hexagonARTMe", false );
        refinementExample( CornerTableAdaptiveRefinement::REDGREEN, "hexagonRedGreen", false );
        refinementExample( CornerTableAdaptiveRefinement::INCREMENTAL, "hexagonIncremental", false );
    }

    //Circle example with three different patterns.
    {
        refinementExample( CornerTableAdaptiveRefinement::ARTMe, "sphereARTMe", true );
        refinementExample( CornerTableAdaptiveRefinement::REDGREEN, "sphereRedGreen", true );
        refinementExample( CornerTableAdaptiveRefinement::INCREMENTAL, "sphereIncremental", true );
    }


    return 0;
}
