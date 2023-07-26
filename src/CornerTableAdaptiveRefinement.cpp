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


#include "CornerTableAdaptiveRefinement.h"
#include <ctime>
#include <cstdio>
#include <cstdlib>
#include <fstream>



void CornerTableAdaptiveRefinement::computeMiddlePoint( const CornerTable* cornerTable, const CornerType c, double* attributes )
{
    //Get the index of two vertices of the edge.
    const CornerType p1 = cornerTable->cornerToVertexIndex( cornerTable->cornerNext( c ) );
    const CornerType p2 = cornerTable->cornerToVertexIndex( cornerTable->cornerPrevious( c ) );

    //Get the vector with geometry information.
    const double* attributesOfVertices = cornerTable->getAttributes( );

    //Get the number of the coordinates by vertex.
    unsigned int n = cornerTable->getNumberAttributesByVertex( );

    //Compute the middle point.
    for (unsigned int i = 0; i < n; i++)
    {
        attributes[i] = ( attributesOfVertices[n * p1 + i] + attributesOfVertices[n * p2 + i] ) / 2.0;
    }
}



CornerTableAdaptiveRefinement::CornerTableAdaptiveRefinement( CornerTable* surface )
{
    //Save a pointer to corner table.
    _surface = surface;

    //Initialize the auxiliar vector.
    _coordinates = 0;

    //Initialize the vector to store the refinement levels.
    _triangleLevels.resize( _surface->getNumTriangles( ), 0 );
    _vertexLevels.resize( _surface->getNumberVertices( ), 0 );
}



CornerTableAdaptiveRefinement::CornerTableAdaptiveRefinement( const char* path )
{
    std::ifstream in( path );
    if (in.fail( ))
    {
        printf( "Error opening the file %s\n", path );
        return;
    }

    CornerType numVer, numTri, numCoord;
    in >> numVer >> numTri >> numCoord;

    double *coord = new double[numCoord * numVer];
    CornerType* trian = new CornerType[3 * numTri];

    _coordinates = 0;

    //Initialize the vectors.
    _triangleLevels.resize( numTri, 0 );
    _vertexLevels.resize( numVer, 0 );

    for (CornerType i = 0; i < numVer; i++)
    {
        for (unsigned int c = 0; c < numCoord; c++)
        {
            in >> coord[2 * i + c];
        }
        int l;
        in >> l;
        _vertexLevels[i] = l;
    }
    for (CornerType i = 0; i < numTri; i++)
    {
        in >> trian[3 * i + 0] >> trian[3 * i + 1] >> trian[3 * i + 2];
        int l, s;
        in >> l >> s;
        setTriangleLevels( &_triangleLevels[i], l, s );
    }

    _surface = new CornerTable( trian, coord, numTri, numVer, numCoord );
    delete []trian;
    delete[]coord;
}



void CornerTableAdaptiveRefinement::writeMesh( const char* path )
{
    std::ofstream out( path );
    if (out.fail( ))
    {
        printf( "Error opening the file %s\n", path );
        return;
    }

    CornerType numPoints = _surface->getNumberVertices( );
    CornerType numberTriangles = _surface->getNumTriangles( );
    unsigned int numAttrib = _surface->getNumberAttributesByVertex( );

    //Initialize the header.
    out << numPoints << " " << numberTriangles << " " << numAttrib << std::endl;

    double* coords = _surface->getAttributes( );
    const CornerType* triangles = _surface->getTriangleList( );

    //Write the coordinates.
    for (CornerType i = 0; i < numPoints; i++)
    {
        for (unsigned int c = 0; c < numAttrib; c++)
        {
            out << coords[numAttrib * i + c] << " ";
        }
        out << ( int ) _vertexLevels[i] << std::endl;
    }

    for (CornerType i = 0; i < numberTriangles; i++)
    {
        out << triangles[3 * i + 0] << " " << triangles[3 * i + 1] << " " << triangles[3 * i + 2] << " ";
        out << ( int ) getTriangleLevel( i ) << " " << ( int ) getTriangleSubLevel( i ) << std::endl;
    }
}



void CornerTableAdaptiveRefinement::writeMeshOFFFIle( const char* path )
{
    std::ofstream out( path );
    if (out.fail( ))
    {
        printf( "Error opening the file %s\n", path );
        return;
    }

    CornerType numPoints = _surface->getNumberVertices( );
    CornerType numberTriangles = _surface->getNumTriangles( );
    unsigned int numAttrib = _surface->getNumberAttributesByVertex( );

    //Write the header.
    out << "OFF\n" << numPoints << " " << numberTriangles << " 0" << std::endl;

    const double* coords = _surface->getAttributes( );
    const CornerType* triangles = _surface->getTriangleList( );

    //Write the coordinates.
    for (CornerType i = 0; i < numPoints; i++)
    {
        for (unsigned int c = 0; c < numAttrib; c++)
        {
            out << coords[numAttrib * i + c] << " ";
        }
        out << std::endl;
    }

    for (CornerType i = 0; i < numberTriangles; i++)
    {
        out << " 3 " << triangles[3 * i + 0] << " " << triangles[3 * i + 1] << " " << triangles[3 * i + 2] << std::endl;
    }
}



CornerTableAdaptiveRefinement::~CornerTableAdaptiveRefinement( )
{
    //Delete auxiliar vector.
    delete _coordinates;

    //Delete queues.
    while (!_undos.empty( ))
    {
        UndoRedoRefinement *undo = _undos.top( );
        delete undo;
        _undos.pop( );
    }

    while (!_redos.empty( ))
    {
        UndoRedoRefinement *redo = _redos.top( );
        delete redo;
        _redos.pop( );
    }
}



CornerTable* CornerTableAdaptiveRefinement::getSurface( )
{
    return _surface;
}



LevelType CornerTableAdaptiveRefinement::getVertexLevel( const CornerType v )
{
    return _vertexLevels[v];
}



void CornerTableAdaptiveRefinement::runSplit( std::vector<bool>& processedTriangles,
                                              std::vector<bool>& visitedTriangles,
                                              std::queue<NodeOfSearch>& cornersForSplit,
                                              const CornerType cornerOperation,
                                              UndoRedoRefinement* undoRedo,
                                              const bool undo,
                                              void (*computeNewVertex )( const CornerTable* cornerTable,
                                                                         const CornerType corner, double *attributes ) )
{
    //Get the opposite corner to the corner that will be used to apply the Edge
    //Split operation.
    CornerType opposite = _surface->cornerOpposite( cornerOperation );

    //Get the number of triangles on current mesh.
    const CornerType numTriangles = _surface->getNumTriangles( );

    //Mark the faces as processed.
    processedTriangles[cornerOperation / 3] = true;
    processedTriangles[numTriangles] = true;

    if (undo)
    {
        CornerType op;
        CornerType c = _surface->cornerNext( cornerOperation );
        setOperation( &op, c, UndoRedoRefinement::WELD );
        undoRedo->addUndo( op );
    }

    //Allocate the vector to store the new coordinates.
    double* coordinates = new double[_surface->getNumberAttributesByVertex( )];

    //Compute new coordinates.
    computeNewVertex( _surface, cornerOperation, coordinates );

    //Perform the edge split operation.
    _surface->edgeSplit( cornerOperation, coordinates );

    //Free the memory.
    delete [] coordinates;

    //New pair of vertex and edge to be queued.
    NodeOfSearch v;
    v._vertex = _surface->cornerNext( cornerOperation );
    v._edge = CornerTable::BORDER_CORNER;

    //If there is an opposite triangle, update the status of the new two opposite
    //triangles.
    if (opposite != CornerTable::BORDER_CORNER)
    {
        //If the triangle was already visited, mark them as processed.
        if (visitedTriangles[opposite / 3])
        {
            processedTriangles[opposite / 3] = true;
            processedTriangles[numTriangles + 1] = true;
        }
        else
        {
            //Set edge to be flipped.
            v._edge = _surface->cornerPrevious( opposite );
        }

        //Mark faces as visited.
        visitedTriangles[opposite / 3] = true;
        visitedTriangles[numTriangles + 1] = true;
    }

    //Add the new pair on queued.
    cornersForSplit.push( v );
}



void CornerTableAdaptiveRefinement::refineTriangles( const CornerType t,
                                                     std::vector<bool>& processedTriangles,
                                                     std::vector<bool>& visitedTriangles,
                                                     std::queue<NodeOfSearch>& cornersForSplit,
                                                     UndoRedoRefinement* undoRedo, const bool undo,
                                                     void (*computeNewVertex )( const CornerTable* cornerTable,
                                                                                const CornerType corner, double *attributes ) )
{
    //Get the initial corner.
    const CornerType initialCorner = 3 * t;

    //Obtain the number of triangles on mesh.
    const CornerType numberTriangles = _surface->getNumTriangles( );

    //Perform the firts split to begin the refinement.
    runSplit( processedTriangles, visitedTriangles, cornersForSplit, initialCorner, undoRedo, undo, computeNewVertex );

    //Mark the triangles as visited.
    visitedTriangles[initialCorner / 3] = true;
    visitedTriangles[numberTriangles] = true;

    //In this initial step mark the triangles as not processed.
    processedTriangles[initialCorner / 3] = false;
    processedTriangles[numberTriangles] = false;

    //Save the edge on the initial triangle to be flipped at the final.
    CornerType cornerFliped = _surface->cornerPrevious( initialCorner );

    NodeOfSearch v;

    //Traverse the mesh refining the triangles.
    while (!cornersForSplit.empty( ))
    {
        //Get the first element on queue.
        v = cornersForSplit.front( );

        //Get the vertex to refine its star.
        CornerType currentCorner = v._vertex;

        //Define the first corner as the current corner.
        CornerType firstCorner = currentCorner;

        //Define the corner to apply the edge split operation.
        CornerType cornerOperation = currentCorner;

        //Refine the star of the vertex.
        cornersForSplit.pop( );
        do
        {
            //Refine just not processed triangles.
            if (!processedTriangles[cornerOperation / 3])
            {
                runSplit( processedTriangles, visitedTriangles, cornersForSplit, cornerOperation, undoRedo, undo, computeNewVertex );
            }

            currentCorner = _surface->cornerRight( cornerOperation );
            cornerOperation = _surface->cornerNext( currentCorner );
        }
        while (cornerOperation != firstCorner && currentCorner != CornerTable::BORDER_CORNER);

        //If the search complete the turn, stop the spinning and perform an Edge
        //Flip operation.
        if (currentCorner != CornerTable::BORDER_CORNER)
        {
            tryFlipEdge( v._edge, undo, undoRedo );
            continue;
        }


        //If stop on a border, spin on the another direction.
        cornerOperation = firstCorner;
        currentCorner = firstCorner;

        do
        {
            //Refine just not processed triangles.
            if (!processedTriangles[cornerOperation / 3])
            {
                runSplit( processedTriangles, visitedTriangles, cornersForSplit, cornerOperation, undoRedo, undo, computeNewVertex );
            }
            currentCorner = _surface->cornerLeft( cornerOperation );
            cornerOperation = _surface->cornerPrevious( currentCorner );
        }
        while (currentCorner != CornerTable::BORDER_CORNER);

        tryFlipEdge( v._edge, undo, undoRedo );
    }

    //Perform the edge flip operation on the initial triangle.
    tryFlipEdge( cornerFliped, undo, undoRedo );
}



void CornerTableAdaptiveRefinement::meshGlobalRefinement( const bool undo, void(*computeNewVertex )( const CornerTable*,
                                                                                                     const CornerType, double* ) )
{
    //Obtain the number of triangles on mesh.
    const CornerType numberTriangles = _surface->getNumTriangles( );

    //Allocate vectors to manage face status.
    std::vector< bool > visitedFaces( 4 * numberTriangles ), processedFaces( 4 * numberTriangles );

    //Allocate queue to traverse on mesh
    std::queue<NodeOfSearch> cornersForSplit;

    //Create an undo object to store the undo operations.
    UndoRedoRefinement *undoRedo = new UndoRedoRefinement( _surface );

    //Refine faces.
    for (CornerType i = 0; i < numberTriangles; i++)
    {
        if (!processedFaces[i])
        {
            refineTriangles( i, processedFaces, visitedFaces, cornersForSplit, undoRedo, undo, computeNewVertex );
        }
    }

    updateUndos( undoRedo );
}



void CornerTableAdaptiveRefinement::meshLocalRefinement( const std::set<CornerType>& trianglesForRefine,
                                                         const bool undo,
                                                         void(*computeNewVertex )( const CornerTable*, const CornerType, double* ) )
{
    //Obtain the number of triangles on mesh.
    const CornerType numberTriangles = _surface->getNumTriangles( );

    //Allocate vectors to manage face status.
    std::vector< bool > visitedTriangles( 4 * numberTriangles, true ), processedTriangles( 4 * numberTriangles );

    //Mark faces to refine as visited.
    for (std::set<CornerType>::const_iterator i = trianglesForRefine.begin( ); i != trianglesForRefine.end( ); i++)
    {
        visitedTriangles[*i] = false;
    }

    //Allocate queue to traverse on mesh
    std::queue<NodeOfSearch> cornersForSplit;

    //Create an undo object to store the undo operations.
    UndoRedoRefinement *undoRedo = new UndoRedoRefinement( _surface );

    //Refine faces.
    for (std::set<CornerType>::const_iterator i = trianglesForRefine.begin( ); i != trianglesForRefine.end( ); i++)
    {
        if (!processedTriangles[*i])
        {
            refineTriangles( *i, processedTriangles, visitedTriangles, cornersForSplit, undoRedo, undo, computeNewVertex );
        }
    }

    updateUndos( undoRedo );
}



void CornerTableAdaptiveRefinement::meshAdaptiveRefinement( const std::vector<CornerType>& trianglesForRefine,
                                                            const bool undo, const Triangulation triangulation,
                                                            const unsigned int ratio,
                                                            void (*computeNewVertex )( const CornerTable* cornerTable,
                                                                                       const CornerType corner, double *attributes ) )
{
    //Set of triangles selected to be refined.
    std::set<Triangle> triangles;

    //Mount the set of faces organized by levels.
    for (CornerType i = 0; i < trianglesForRefine.size( ); i++)
    {
        Triangle f;
        f._index = trianglesForRefine[i];
        f._level = getTriangleLevel( f._index );
        triangles.insert( f );
    }

    //Perform the suitable preprocessing step.
    if (triangulation == ARTMe)
    {
        preprocessingARTMe( triangles );
    }
    else if (triangulation == REDGREEN)
    {
        preprocessingRedGreen( triangles );
    }
    else if (triangulation == INCREMENTAL)
    {
        preprocessingIncremental( triangles, ratio );
    }
    else
    {
        printf( "Invalid kind of triangulation\n" );
        return;
    }

    //Refine the triangles.
    UndoRedoRefinement *undoRedo = makeRefinement( triangles, undo, computeNewVertex );

    //Update the undo stack.
    updateUndos( undoRedo );
}



void CornerTableAdaptiveRefinement::meshAdaptiveRefinement( const std::set<CornerType>& trianglesForRefine,
                                                            const bool undo, const Triangulation triangulation,
                                                            const unsigned int ratio,
                                                            void (*computeNewVertex )( const CornerTable* cornerTable,
                                                                                       const CornerType corner, double *attributes ) )
{
    //Set of triangles selected to be refined.
    std::set<Triangle> triangles;

    //Mount the set of faces organized by levels.
    for (std::set<CornerType>::iterator it = trianglesForRefine.begin( );
         it != trianglesForRefine.end( ); it++)
    {
        Triangle f;
        f._index = *it;
        f._level = getTriangleLevel( f._index );
        triangles.insert( f );
    }

    //Perform the suitable preprocessing step.
    if (triangulation == ARTMe)
    {
        preprocessingARTMe( triangles );
    }
    else if (triangulation == REDGREEN)
    {
        preprocessingRedGreen( triangles );
    }
    else if (triangulation == INCREMENTAL)
    {
        preprocessingIncremental( triangles, ratio );
    }
    else
    {
        printf( "Invalid kind of triangulation\n" );
        return;
    }

    //Refine the triangles.
    UndoRedoRefinement *undoRedo = makeRefinement( triangles, undo, computeNewVertex );

    //Update the undo stack.
    updateUndos( undoRedo );
}



void CornerTableAdaptiveRefinement::preprocessingARTMe( std::set<Triangle>& trianglesForRefine )
{
    //Vector with  the new set of triangles to be refined.
    std::vector<Triangle> newTriangles;

    //Copy the set of triangles to the new vector.
    newTriangles.insert( newTriangles.end( ), trianglesForRefine.begin( ), trianglesForRefine.end( ) );

    //Add neighbors with lower level and replace the triangles with sublevel 3.
    for (CornerType i = 0; i < newTriangles.size( ); i++)
    {
        //Get the triangle index.
        CornerType t = newTriangles[i]._index;

        //Get the triangle sublevel.
        unsigned char sublevel = getTriangleSubLevel( newTriangles[i]._index );

        if (sublevel != 3)
        {
            //Add neighbor triangles with a lower level to be refined in order to
            //generate a balanced mesh.
            addLowerLevelTriangles( trianglesForRefine, newTriangles, newTriangles[i] );
        }
        else
        {
            //Get the neighbor (sibling) triangle with sublevel 3.
            Triangle f;
            if (addNeighbor3SublevelTriangles( t, f ))
            {
                //Remove the triangle with sublevel 3.
                trianglesForRefine.erase( newTriangles[i] );

                if (trianglesForRefine.find( f ) == trianglesForRefine.end( ))
                {
                    //Add the sibling triangle to be refined.
                    trianglesForRefine.insert( f );

                    //Add the sibling triangle to be preprocessed too.
                    newTriangles.push_back( f );
                }
            }
            else
            {
                printf( "Sibling triangle of sublevel 3 does not found\n" );
                return;
            }
        }
    }
}



bool CornerTableAdaptiveRefinement::addNeighbor3SublevelTriangles( CornerType triangle, Triangle &t )
{
    //Get the triangle level.
    LevelType level = getTriangleLevel( triangle );

    //Visit the neighbor triangles.
    for (CornerType corner = 3 * triangle; corner != 3 * triangle + 3; corner++)
    {
        //Get a corner on the opposite triangle.
        CornerType oppositeTriangle = _surface->cornerOpposite( corner );

        //Verify if that triangle exists.
        if (oppositeTriangle != CornerTable::BORDER_CORNER)
        {
            //Get the triangle index.
            oppositeTriangle /= 3;

            //Get the sublevel of the oposite face.
            unsigned char oppositeSubLevel = getTriangleSubLevel( oppositeTriangle );

            //Do something just if the sublevel different of 3.
            if (oppositeSubLevel != 3)
            {
                //Get the corners of the opposite edge.
                CornerType next = _surface->cornerNext( corner );
                CornerType prev = _surface->cornerPrevious( corner );

                //Get the vertices index of the edge.
                CornerType v1 = _surface->cornerToVertexIndex( next );
                CornerType v2 = _surface->cornerToVertexIndex( prev );

                //Verify if the vertices level is greater than the triangle level.
                //If it is true, this edge divides the current triangle and a
                //triangle with sublevel 0 in the next level.
                if (_vertexLevels[v1] == level + 1 && _vertexLevels[v2] == level + 1)
                {
                    //Get the index of the triangle to be refined.
                    CornerType triangleToRefined = CornerTable::BORDER_CORNER;

                    CornerType t1 = _surface->cornerRight( corner );
                    CornerType t2 = _surface->cornerLeft( corner );

                    if (t1 != CornerTable::BORDER_CORNER)
                    {
                        t1 /= 3;

                        if (getTriangleSubLevel( t1 ) == 1 && getTriangleLevel( t1 ) == level)
                        {
                            triangleToRefined = t1;
                        }
                    }
                    if (t2 != CornerTable::BORDER_CORNER)
                    {
                        t2 /= 3;
                        if (getTriangleSubLevel( t2 ) == 2 && getTriangleLevel( t2 ) == level)
                        {
                            triangleToRefined = t2;
                        }
                    }
                    if (triangleToRefined == CornerTable::BORDER_CORNER)
                    {
                        printf( "Sibling triangle of sublevel 3 does not found\n" );
                        return false;
                    }

                    //Set informations of the triangle to be refined.
                    t._index = triangleToRefined;
                    t._level = getTriangleLevel( t._index );

                    return true;
                }
            }
        }
    }
    return false;
}



void CornerTableAdaptiveRefinement::addLowerLevelTriangles( std::set<Triangle>& trianglesForRefine,
                                                            std::vector<Triangle>& newTriangles,
                                                            const Triangle& t )
{
    //Get the triangle index.
    CornerType triangle = t._index;

    //Get the triangle level.
    LevelType level = t._level;

    for (CornerType corner = 3 * triangle; corner != 3 * triangle + 3; corner++)
    {
        //Get a corner in the opposite triangle.
        CornerType oppositeTriangle = _surface->cornerOpposite( corner );

        //Verify is this triangle exists.
        if (oppositeTriangle != CornerTable::BORDER_CORNER)
        {
            //Get the triangle index..
            oppositeTriangle /= 3;

            //Get the level of the opposite triangle.
            LevelType oppositeLevel = getTriangleLevel( oppositeTriangle );

            //Verify if the opposite triangle has a lower refinement level. If
            //it is true, this triangle is add to be refined.
            if (oppositeLevel < level)
            {
                //Mount the triangle informations.
                Triangle f;
                f._index = oppositeTriangle;
                f._level = getTriangleLevel( f._index );

                //Verify is this triangle has been already added.
                if (trianglesForRefine.find( f ) == trianglesForRefine.end( ))
                {
                    if (getTriangleSubLevel( oppositeTriangle ) != 3)
                    {
                        //Add to be refined just if the triangle have a refinement
                        //level different of 3.
                        trianglesForRefine.insert( f );
                    }

                    //Add new triangle to be preprocessed.
                    newTriangles.push_back( f );
                }
            }
        }
    }
}



void CornerTableAdaptiveRefinement::addSiblingTriangles( std::set<Triangle>& trianglesForRefine,
                                                         std::vector<Triangle>& newTriangles,
                                                         std::vector<bool>& visitedTriangles,
                                                         const Triangle& t )
{
    //Get the triangle index.
    CornerType triangle = t._index;

    //Get the triangle level..
    LevelType level = t._level;

    //Get the triangle sublevel.
    LevelType sublevel = getTriangleSubLevel( triangle );

    //Just do something if the current triangle has sublevel 1 or 2. If it has
    //a sublevel 1, the sibling triangle with sublevel 2 will be added to be
    //refinement and vice-versa.
    if (sublevel == 1 || sublevel == 2)
    {
        //Search by the vertex with greatest level on face.
        for (CornerType corner = 3 * triangle; corner != 3 * triangle + 3; corner++)
        {
            CornerType v = _surface->cornerToVertexIndex( corner );
            if (level + 1 == getVertexLevel( v ))
            {
                //Get the sibling triangle to be tested.
                CornerType oppositeTriangle = sublevel == 1 ? _surface->cornerRight( corner ) : _surface->cornerLeft( corner );
                LevelType siblingTriangleSublevel = sublevel == 1 ? 2 : 1;

                //Mount the triangle information.
                Triangle f;
                f._index = oppositeTriangle / 3;
                f._level = getTriangleLevel( f._index );

                //Verify if the sibling triangle was found.
                if (getTriangleSubLevel( f._index ) == siblingTriangleSublevel && getTriangleLevel( f._index ) == level &&
                    trianglesForRefine.find( f ) == trianglesForRefine.end( ))
                {
                    //Insert the new triangle to be refined and to be preprocessed.
                    trianglesForRefine.insert( f );
                    newTriangles.push_back( f );
                }

                //Do the verification to the opposite triangle.
                oppositeTriangle = _surface->cornerOpposite( corner );
                if (oppositeTriangle != CornerTable::BORDER_CORNER)
                {
                    f._index = oppositeTriangle / 3;
                    f._level = getTriangleLevel( f._index );

                    //Verify if the sibling triangle was found.
                    if (getTriangleSubLevel( f._index ) == siblingTriangleSublevel && getTriangleLevel( f._index ) == level
                        && trianglesForRefine.find( f ) == trianglesForRefine.end( ))
                    {
                        trianglesForRefine.insert( f );
                        newTriangles.push_back( f );
                    }
                    else if (getTriangleSubLevel( f._index ) == 0 && trianglesForRefine.find( f ) == trianglesForRefine.end( ))
                    {
                        //Verify with the triangle will has two or more cracks.
                        if (visitedTriangles[f._index])
                        {
                            trianglesForRefine.insert( f );
                            newTriangles.push_back( f );
                        }
                        visitedTriangles[f._index] = true;
                    }
                }
                break;
            }
        }
    }
}



void CornerTableAdaptiveRefinement::addOpositeTriangles( std::set<Triangle>& trianglesForRefine,
                                                         std::vector<Triangle>& newTriangles,
                                                         std::vector<bool>& visitedTriangles,
                                                         const Triangle& t )
{
    //Get the triangle index.
    CornerType triangle = t._index;

    //Get the triangle level.
    LevelType level = t._level;

    //Get the triangle sublevel.
    LevelType sublevel = getTriangleSubLevel( triangle );

    for (CornerType corner = 3 * triangle; corner != 3 * triangle + 3; corner++)
    {
        //Get a corner on the opposite triangle.
        CornerType oppositeTriangle = _surface->cornerOpposite( corner );

        //Verify if the triangle exists.
        if (oppositeTriangle != CornerTable::BORDER_CORNER)
        {
            //Get the index of the opposite triangle.
            oppositeTriangle /= 3;

            //Get the sublevel of the opposite triangle.
            LevelType oppositeSubLevel = getTriangleSubLevel( oppositeTriangle );

            //Mount the triangle information.
            Triangle f;
            f._index = oppositeTriangle;
            f._level = getTriangleLevel( f._index );

            //Verify if the opposite triangle has two or more cracks.
            if (sublevel == 0 && oppositeSubLevel == 0 &&
                trianglesForRefine.find( f ) == trianglesForRefine.end( ))
            {
                if (visitedTriangles[oppositeTriangle])
                {
                    trianglesForRefine.insert( f );
                    newTriangles.push_back( f );
                }
                visitedTriangles[oppositeTriangle] = true;
            }
                //In case of the opposite case has a sublevel different of 0, a
                //triangle with sublevel 3 will be generated, thus add this face
                //to be refined.
            else if (oppositeSubLevel != 0)
            {
                if (getTriangleLevel( oppositeTriangle ) == level && trianglesForRefine.find( f ) == trianglesForRefine.end( ))
                {
                    trianglesForRefine.insert( f );
                    newTriangles.push_back( f );
                }
            }
        }
    }
}



void CornerTableAdaptiveRefinement::preprocessingRedGreen( std::set<Triangle>& trianglesForRefine )
{
    //New faces to be refined.
    std::vector<Triangle> newTriangles;

    //A vector to speed up the preprocessing process.
    std::vector<bool> visitedTriangles( _surface->getNumTriangles( ), false );

    //Copy the selected triangles to the new vector.
    newTriangles.insert( newTriangles.end( ), trianglesForRefine.begin( ), trianglesForRefine.end( ) );

    //Trata faces vizinhas com niveis diferentes onde uma foi marcada para 
    //subdivisao e outra nao. No caso de a nao marcada ter um nivel de subdivisao
    //menor, ambas devem ser marcadas para subdivisao.

    //Verify if there are trianhhles with lower refinement level and triangles that
    //will generate a triangle with sublevel 3.
    for (CornerType i = 0; i < newTriangles.size( ); i++)
    {
        //Add triangles with a lower level.
        addLowerLevelTriangles( trianglesForRefine, newTriangles, newTriangles[i] );

        //Add sibling triangles.
        addSiblingTriangles( trianglesForRefine, newTriangles, visitedTriangles, newTriangles[i] );

        //Add opposite triangles.
        addOpositeTriangles( trianglesForRefine, newTriangles, visitedTriangles, newTriangles[i] );
    }
}



void CornerTableAdaptiveRefinement::preprocessingIncremental( std::set<Triangle>& trianglesForSubdivide, unsigned int radius )
{
    //Corners of the vertices to do the expansion step.
    std::set<CornerType> cornersForExpansion, aux;

    //Get the corners of the vertices to do the expansion step.
    for (std::set<Triangle>::iterator i = trianglesForSubdivide.begin( ); i != trianglesForSubdivide.end( ); i++)
    {
        for (CornerType c = 3 * i->_index; c < 3 * i->_index + 3; c++)
        {
            CornerType v = _surface->cornerToVertexIndex( c );

            if (aux.find( v ) == aux.end( ))
            {
                //Adiciona corner na lista
                cornersForExpansion.insert( c );
                aux.insert( v );
            }
        }
    }
    aux.clear( );

    //A set to get the corner of that vertex on the expanded region.
    std::set<CornerType> vertexsForExpansion;

    //Compute the expansion of radius r.
    getRRingFromVertices( cornersForExpansion, radius, vertexsForExpansion );

    //Compute that triangles formed by the vertices on the expanded region..
    computeTrianglesFormedByVertexs( vertexsForExpansion, trianglesForSubdivide );

    //Perform the preprocessing ARTMe. It is necessary just on hybrid cases.
    preprocessingARTMe( trianglesForSubdivide );
}



void CornerTableAdaptiveRefinement::computeTrianglesFormedByVertexs( const std::set<CornerType>& vertices,
                                                                     std::set<Triangle>& triangles )
{
    //Get the number of points.
    CornerType numberPoints = _surface->getNumberVertices( );

    //Vector to mark that points on the expanded region.
    bool *verticesForRefinement = new bool[ numberPoints ];

    //Initialize the vector.
    memset( verticesForRefinement, 0, numberPoints * sizeof ( bool ) );

    //Mark that vertices on the selected are.
    for (std::set<CornerType>::iterator it = vertices.begin( ); it != vertices.end( ); ++it)
    {
        verticesForRefinement[*it] = true;
    }

    //Search by the triangles on the expanded are to be refined..
    for (std::set<CornerType>::iterator it = vertices.begin( ); it != vertices.end( ); ++it)
    {
        //Obtem corner do vertice corrente.
        const CornerType corner = _surface->vertexToCornerIndex( *it );

        //Traverse the star of the vertex.
        CornerType currentCorner = _surface->cornerNext( corner );
        const CornerType firstCorner = currentCorner;
        do
        {
            //Get the the vertex of the current edge opposite to the current
            //corner.
            CornerType vertexA = _surface->cornerToVertexIndex( currentCorner );
            CornerType vertexB = _surface->cornerToVertexIndex( _surface->cornerNext( currentCorner ) );

            //If this two vertices are in the expanded set, there is a new triangle
            //to be refined.
            if (verticesForRefinement[ vertexA ] && verticesForRefinement[vertexB ])
            {
                //Mount the triangle information.
                Triangle f;
                f._index = currentCorner / 3;
                f._level = getTriangleLevel( f._index );

                //Add this new triangle on the set to be refined.
                triangles.insert( f );
            }
            currentCorner = _surface->cornerRight( currentCorner );
        }
        while (currentCorner != CornerTable::BORDER_CORNER && currentCorner != firstCorner);

        //Stop if complete the spin.
        if (currentCorner != CornerTable::BORDER_CORNER)
            continue;

        //Spin in the opposite orientation.
        currentCorner = _surface->cornerPrevious( corner );
        do
        {
            //Get the the vertex of the current edge opposite to the current
            //corner.
            CornerType vertexA = _surface->cornerToVertexIndex( currentCorner );
            CornerType vertexB = _surface->cornerToVertexIndex( _surface->cornerPrevious( currentCorner ) );

            //If this two vertices are in the expanded set, there is a new triangle
            //to be refined.
            if (verticesForRefinement[vertexA] && verticesForRefinement[vertexB ])
            {
                //Mount the triangle information.
                Triangle f;
                f._index = currentCorner / 3;
                f._level = getTriangleLevel( f._index );

                //Add this new triangle on the set to be refined.
                triangles.insert( f );
            }

            currentCorner = _surface->cornerLeft( currentCorner );
        }
        while (currentCorner != CornerTable::BORDER_CORNER);
    }

    //Free memory.
    delete[] verticesForRefinement;
}



void CornerTableAdaptiveRefinement::getRRingFromVertices( const std::set<CornerType>& corners,
                                                          const unsigned int radius,
                                                          std::set<CornerType>& vertices )
{
    //Get the number of points on the mesh.
    CornerType numberPoints = _surface->getNumberVertices( );

    //Vector to mark a vertex as visited or not.
    bool *vertexVisited = new bool[numberPoints];

    //Store the distance between the region and the vertex.
    unsigned int* mapDistanceVertex = new unsigned int[numberPoints];

    //Pair to store the corner index and the distance.
    std::pair<CornerType, int> cornerDistance;

    //Initialize the distance vector.
    memset( mapDistanceVertex, 255, numberPoints * sizeof ( unsigned int ) );

    //Initialize the visited vector.
    memset( vertexVisited, 0, numberPoints * sizeof ( bool ) );

    //List to expand the region..
    std::queue< std::pair<CornerType, int> > listForSearch;

    for (std::set<CornerType>::const_iterator v = corners.begin( ); v != corners.end( ); v++)
    {
        //Get the initial corner.
        CornerType corner = *v;

        //Define the distance variable.
        cornerDistance.first = corner;
        cornerDistance.second = 0;

        //Add the object on the list.
        listForSearch.push( cornerDistance );

        //Get the vertex index.
        CornerType vertex = _surface->cornerToVertexIndex( corner );

        //Insert the vertex on the list.
        vertices.insert( vertex );

        //Mark the vertex as visited.
        vertexVisited[ vertex ] = true;

        //Mark the distance to 0.
        mapDistanceVertex[vertex] = 0;

        //Do the search from this vertex.
        while (!listForSearch.empty( ))
        {
            //Get the first element on the queue.
            cornerDistance = listForSearch.front( );

            //Remove the first element.
            listForSearch.pop( );

            //Save the corner and current distance.
            CornerType currentCorner = cornerDistance.first;
            CornerType currentDistance = cornerDistance.second;

            //Get the neighborhood.
            const std::vector<CornerType> neigbours = _surface->getCornerNeighbours( currentCorner );
            for (unsigned int j = 0; j < neigbours.size( ); j++)
            {
                //Get the vertex index.
                vertex = _surface->cornerToVertexIndex( neigbours[j] );

                //Verify if the vertex was not visited or the distance is greater.
                if (!vertexVisited[vertex] || ( vertexVisited[vertex] && mapDistanceVertex[vertex] > currentDistance + 1 ))
                {
                    //Mark the vertex as visited
                    vertexVisited[vertex] = true;

                    //Verify if the distance is less than the radius.
                    if (currentDistance + 1 <= radius)
                    {
                        //Put the element on queue.
                        cornerDistance.first = neigbours[j];
                        cornerDistance.second = currentDistance + 1;
                        listForSearch.push( cornerDistance );

                        //Addt the vertex on the list of vertices on the region.
                        vertices.insert( vertex );

                        //Update the distance.
                        mapDistanceVertex[vertex] = currentDistance + 1;
                    }
                }
            }
        }
    }

    //Free memory.
    delete[] vertexVisited;
    delete[] mapDistanceVertex;
}



UndoRedoRefinement * CornerTableAdaptiveRefinement::makeRefinement( const std::set<Triangle>& trianglesForRefine,
                                                                    const bool undo,
                                                                    void (*computeNewVertex )( const CornerTable* cornerTable,
                                                                                               const CornerType corner, double *attributes ) )
{
    //Allocate the auxiliary vector to store the new coordinates.
    delete [] _coordinates;
    _coordinates = new double[_surface->getNumberAttributesByVertex( )];

    //Allocate the structure to store the undo operations.
    UndoRedoRefinement *undoRedo = new UndoRedoRefinement( _surface );

    //Get the number of triangles.
    const CornerType numberTriangles = _surface->getNumTriangles( );

    //Vector to store if a triangle is processed or not.
    std::vector<bool> processedTriangles;
    processedTriangles.resize( 4 * numberTriangles, false );

    //Mark all triangles on surface as processed.
    for (CornerType i = 0; i < numberTriangles; i++)
    {
        processedTriangles[i] = true;
    }

    //Mark that triangles that were selected to be refined as unprocessed.
    for (std::set<Triangle>::iterator it = trianglesForRefine.begin( );
         it != trianglesForRefine.end( ); it++)
    {
        processedTriangles[( *it )._index] = false;
    }

    //Refine the triangle from to lower to highest refinement level.
    for (std::set<Triangle>::iterator it = trianglesForRefine.begin( );
         it != trianglesForRefine.end( ); it++)
    {
        //Get the triangle index.
        const CornerType triangle = ( *it )._index;

        //Get the sublevel of the triangle.
        const LevelType subLevel = getTriangleSubLevel( triangle );

        //Verify if the triangle was already processed.
        if (!processedTriangles[triangle])
        {
            //Verify the sublevel of the triangle.
            switch (subLevel)
            {
                case 0:
                {
                    //For the sublevel 0, begin to refine from any vertex.
                    meshRefineTriangles( processedTriangles, 3 * triangle, undo, undoRedo, computeNewVertex );
                    break;
                }

                case 1:
                case 2:
                {
                    //For the sublevel 1 or 2 begin to refine from the vertex
                    //of the highest level.
                    CornerType initialCorner = CornerTable::BORDER_CORNER;
                    for (CornerType corner = 3 * triangle; corner < 3 * triangle + 3; corner++)
                    {
                        //Obtem o indice do vertice corrente.
                        const CornerType vertice = _surface->cornerToVertexIndex( corner );

                        if (_vertexLevels[vertice] == getTriangleLevel( triangle ) + 1)
                        {
                            initialCorner = corner;
                            break;
                        }
                    }
                    if (initialCorner < 0)
                    {
                        printf( "Vertex with high level was not found\n" );
                    }

                    meshRefineTriangles( processedTriangles, initialCorner, undo, undoRedo, computeNewVertex );
                    break;
                }
                case 3:
                    printf( "Error: trying to refined a triangle with sublevel 3\n" );
                    return undoRedo;
                    break;
                default:
                {
                    printf( "Error: sublevel out of [0,3]\n" );
                    return undoRedo;
                    break;
                }
            }
        }
    }

    //Free memory of auxiliary vector.
    delete [] _coordinates;
    _coordinates = 0;

    //Return the stack of undo.
    return undoRedo;
}



void CornerTableAdaptiveRefinement::meshRefineTriangles( std::vector<bool>& processedTriangles,
                                                         const CornerType initialCorner,
                                                         const bool undo,
                                                         UndoRedoRefinement* currentUndoRedo,
                                                         void (*computeNewVertex )( const CornerTable* cornerTable,
                                                                                    const CornerType corner, double *attributes ) )
{
    //Queue for store the elements on search.
    std::queue<NodeOfSearch> cornersForSplit;

    //Get the triangle sublevel.
    const LevelType subLevelTriangle = getTriangleSubLevel( initialCorner / 3 );

    //If it is the first Edge Split on triangle, save the edge corner to be
    //flipped at the final. This must be happens just for sublevel 0..
    CornerType cornerFliped = CornerTable::BORDER_CORNER;
    if (subLevelTriangle == 0)
    {
        cornerFliped = _surface->cornerPrevious( initialCorner );
    }

    //Initialize the first by an Edge Split Operation.
    runSplitAdaptiveRefinement( processedTriangles, cornersForSplit,
                                initialCorner, currentUndoRedo, undo,
                                computeNewVertex );

    //Verify if its necessary to flip the edge on initial triangle. This just
    //can happen when the algorithm begins in a triangle with sublevel 1 or 2.
    CornerType flipCornerOnTriangle = CornerTable::BORDER_CORNER;

    if (subLevelTriangle > 0)
    {
        if (subLevelTriangle == 1)
        {
            flipCornerOnTriangle = _surface->cornerNext( initialCorner );
        }
        else if (subLevelTriangle == 2)
        {
            flipCornerOnTriangle = _surface->cornerLeft( initialCorner );
            if (flipCornerOnTriangle != CornerTable::BORDER_CORNER)
            {
                flipCornerOnTriangle = _surface->cornerNext( flipCornerOnTriangle );
            }
        }
        else
        {
            printf( "Error processing flip on the initial triangle\n" );
            exit( 0 );
        }
    }

    //Try to execute the flip operation.
    tryFlipEdgeAdaptiveRefinement( flipCornerOnTriangle, undo, currentUndoRedo );

    NodeOfSearch v;

    //Traverse the mesh performing the refinement.
    while (!cornersForSplit.empty( ))
    {
        //Get the first element on queue.
        v = cornersForSplit.front( );

        //Get the corner of the vertex.
        CornerType currentCorner = v._vertex;

        //Define the corner as the first.
        CornerType firstCorner = currentCorner;

        //Define the current corner as the corner to execute the operation.
        CornerType cornerOperation = currentCorner;

        //Remove the fist element of the queue.
        cornersForSplit.pop( );
        do
        {
            //Get the sublevel of the current triangle.
            const LevelType subLevel = getTriangleSubLevel( cornerOperation / 3 );

            //Verify if the triangle needs of more refinement.
            if (!processedTriangles[cornerOperation / 3] && ( subLevel == 1 || subLevel == 2 ))
            {
                runSplitAdaptiveRefinement( processedTriangles, cornersForSplit,
                                            cornerOperation, currentUndoRedo, undo,
                                            computeNewVertex );
            }
            currentCorner = _surface->cornerRight( cornerOperation );
            cornerOperation = _surface->cornerNext( currentCorner );
        }
        while (cornerOperation != firstCorner && currentCorner != CornerTable::BORDER_CORNER);

        //If complete the spin try to flip the edge.
        if (currentCorner != CornerTable::BORDER_CORNER)
        {
            tryFlipEdgeAdaptiveRefinement( v._edge, undo, currentUndoRedo );
            continue;
        }

        //If stop on border spin in another direction.
        cornerOperation = firstCorner;
        currentCorner = firstCorner;
        do
        {
            //Get the sublevel of the current triangle.
            const LevelType subLevel = getTriangleSubLevel( cornerOperation / 3 );

            //Verify if the triangle needs of more refinement.
            if (!processedTriangles[cornerOperation / 3] && ( subLevel == 1 || subLevel == 2 ))
            {
                runSplitAdaptiveRefinement( processedTriangles, cornersForSplit,
                                            cornerOperation, currentUndoRedo, undo, computeNewVertex );
            }
            currentCorner = _surface->cornerLeft( cornerOperation );
            cornerOperation = _surface->cornerPrevious( currentCorner );
        }
        while (currentCorner != CornerTable::BORDER_CORNER);

        //Try to flip the edge after processing the vertex.
        tryFlipEdgeAdaptiveRefinement( v._edge, undo, currentUndoRedo );
    }

    //Try to flip the first edge when the initial triangle has sublevel 0.
    tryFlipEdgeAdaptiveRefinement( cornerFliped, undo, currentUndoRedo );
}



void CornerTableAdaptiveRefinement::processTriangleLevel( CornerType t1, CornerType t2 )
{
    //Get the sublevel of the original triangle.
    const LevelType sublevel = getTriangleSubLevel( t1 );

    //Set the refinement level information of the new triangle equal to the old
    //triangle.
    _triangleLevels.push_back( _triangleLevels[t1] );

    switch (sublevel)
    {
        case 0:
        {
            //For sublevel 0, the old face will have sublevel 1 and the new face
            //sublevel 2.
            incrementSubLevel( &_triangleLevels[t1] );
            incrementSubLevel( &_triangleLevels[t2], 2 );
        }
            break;
        case 1:
        {
            //For sublevel 1, the old face will have the sublevel 3 and the new
            //face will have sublevel 0 in the next refinement level.
            incrementSubLevel( &_triangleLevels[t1], 2 );
            incrementSubLevel( &_triangleLevels[t2], 3 );
        }
            break;
        case 2:
        {
            //For sublevel 2, the oold face will have sublevel 0 in the next
            //refinement level and the new face will have sublevel 3.
            incrementSubLevel( &_triangleLevels[t1], 2 );
            incrementSubLevel( &_triangleLevels[t2] );
        }
            break;
        default:
        {
            printf( "Error: trying to process an invalid sublevel in Edge Split Operation\n" );
        }
    }
}



void CornerTableAdaptiveRefinement::processTriangles( std::vector<bool>& processedTriangles,
                                                      CornerType t1, CornerType t2 )
{
    //Get the sublevel of the triangle.
    unsigned char subLevel = getTriangleSubLevel( t1 );

    //If the triangle was visited, mark it as processed.
    if (subLevel < 1 || subLevel > 2)
    {
        processedTriangles[t1] = true;
    }
    processedTriangles[t2] = processedTriangles[t1];
}



void CornerTableAdaptiveRefinement::runSplitAdaptiveRefinement( std::vector<bool>& processedTriangles,
                                                                std::queue<NodeOfSearch>& cornersForSplit,
                                                                const CornerType cornerOperation,
                                                                UndoRedoRefinement* undoRedo,
                                                                const bool undo,
                                                                void (*computeNewVertex )( const CornerTable* cornerTable,
                                                                                           const CornerType corner, double *attributes ) )
{
    //Get the opposite corner.
    const CornerType opposite = _surface->cornerOpposite( cornerOperation );

    //Get the number of triangles on current mesh.
    const CornerType numTriangles = _surface->getNumTriangles( );

    //Set the level of the new vertex as the triangle level plus one.
    _vertexLevels.push_back( getTriangleLevel( cornerOperation / 3 ) + 1 );

    //Update ethe level of the triangles.
    processTriangleLevel( cornerOperation / 3, numTriangles );

    //Update the status of processed triangles.
    processTriangles( processedTriangles, cornerOperation / 3, numTriangles );

    if (undo)
    {
        //Add undo operation.
        CornerType op;
        CornerType c = _surface->cornerNext( cornerOperation );
        setOperation( &op, c, UndoRedoRefinement::WELD );
        undoRedo->addUndo( op );
    }

    //Verify if it is necessary to flip the edge on the opposite triangle.
    CornerType flipCornerOnOppositeTriangle = CornerTable::BORDER_CORNER;

    if (opposite != CornerTable::BORDER_CORNER)
    {
        //Get the triangle sublevel.
        LevelType subLevelOppositetriangle = getTriangleSubLevel( opposite / 3 );
        if (subLevelOppositetriangle == 1)
        {
            flipCornerOnOppositeTriangle = _surface->cornerOpposite( cornerOperation );
            if (flipCornerOnOppositeTriangle != CornerTable::BORDER_CORNER)
            {
                flipCornerOnOppositeTriangle = _surface->cornerNext( flipCornerOnOppositeTriangle );
            }
        }
    }

    //Compute new coordinates.
    computeNewVertex( _surface, cornerOperation, _coordinates );

    //Perform the Edge Split operation.
    _surface->edgeSplit( cornerOperation, _coordinates );

    //Mount the information to enqueue.
    NodeOfSearch v;
    v._vertex = _surface->cornerNext( cornerOperation );
    v._edge = CornerTable::BORDER_CORNER;

    //Verify if there is an opposite triangle.
    if (opposite != CornerTable::BORDER_CORNER)
    {
        LevelType subLevelOppositeTriangle = getTriangleSubLevel( opposite / 3 );

        //If the sublevel is 0, save the edge to be flipped after processing the
        //vertex.
        if (subLevelOppositeTriangle == 0)
        {
            if (!processedTriangles[opposite / 3])
            {
                v._edge = _surface->cornerPrevious( opposite );
            }
        }

        //Update ethe level of the triangles.
        processTriangleLevel( opposite / 3, numTriangles + 1 );

        //Update the status of processed triangles.
        processTriangles( processedTriangles, opposite / 3, numTriangles + 1 );

        ///Try to flip the edge on the opposite triangle.
        if (subLevelOppositeTriangle == 2)
        {
            flipCornerOnOppositeTriangle = _surface->cornerOpposite( cornerOperation );
            if (flipCornerOnOppositeTriangle != CornerTable::BORDER_CORNER)
            {
                flipCornerOnOppositeTriangle = _surface->cornerPrevious( flipCornerOnOppositeTriangle );
            }
        }
        else if (subLevelOppositeTriangle > 2)
        {
            printf( "Error: invalid sublevel processing the Edge Flip operation\n" );
            return;
        }
        tryFlipEdgeAdaptiveRefinement( flipCornerOnOppositeTriangle, undo, undoRedo );
    }

    //Add element on queue.
    cornersForSplit.push( v );
}



void CornerTableAdaptiveRefinement::tryFlipEdge( const CornerType corner, bool undo, UndoRedoRefinement * undoRedo )
{
    if (corner == CornerTable::BORDER_CORNER || corner == CornerTable::BORDER_CORNER)
    {
        return;
    }

    //Get the opposite corner.
    CornerType opposite = _surface->cornerOpposite( corner );

    //Verify is there is an opposite triangle.
    if (opposite != CornerTable::BORDER_CORNER)
    {
        //Flip the edge.
        _surface->edgeFlip( corner );

        //Add inverse operations.
        if (undo)
        {
            CornerType op;
            CornerType c = _surface->cornerNext( corner );
            setOperation( &op, c, UndoRedoRefinement::UNFLIP );
            undoRedo->addUndo( op );
        }
    }
}



void CornerTableAdaptiveRefinement::tryFlipEdgeAdaptiveRefinement( const CornerType corner,
                                                                   bool undo,
                                                                   UndoRedoRefinement * undoRedo )
{
    if (corner == CornerTable::BORDER_CORNER || corner == CornerTable::BORDER_CORNER)
    {
        return;
    }

    //Get the opposite corner.
    CornerType opposite = _surface->cornerOpposite( corner );

    //Verify is there is an opposite triangle.
    if (opposite != CornerTable::BORDER_CORNER)
    {
        //Verify if it is possible to flip.
        if (getTriangleSubLevel( opposite / 3 ) == 3 && getTriangleSubLevel( corner / 3 ) == 3)
        {
            //Flip the edge.
            _surface->edgeFlip( corner );

            //Update de level of the triangles.
            incrementSubLevel( &_triangleLevels[corner / 3] );
            incrementSubLevel( &_triangleLevels[opposite / 3] );

            //Add inverse operations.
            if (undo)
            {
                CornerType op;
                CornerType c = _surface->cornerNext( corner );
                setOperation( &op, c, UndoRedoRefinement::UNFLIP );
                undoRedo->addUndo( op );
            }
        }
    }
}



void CornerTableAdaptiveRefinement::updateUndos( const UndoRedoRefinement * undo )
{
    if (undo->isEmpty( ))
    {
        delete undo;
        return;
    }
    
    //Add the undo operations on the stack of undo.
    _undos.push( ( UndoRedoRefinement* ) undo );

    //Free the redo stack.
    while (!_redos.empty( ))
    {
        UndoRedoRefinement* redo = _redos.top( );
        delete redo;
        _redos.pop( );
    }
}



void CornerTableAdaptiveRefinement::makeRedo( )
{
    //Do nothing is the stack is empty.
    if (_redos.empty( ))
    {
        return;
    }

    //Get the top of the stack.
    UndoRedoRefinement* redo = _redos.top( );

    //Remove the top of the stack.
    _redos.pop( );

    //Perform the redo operations.
    while (redo->makeRedo( _triangleLevels, _vertexLevels ));

    //Add the redo operations on stack of undo.
    _undos.push( redo );
}



void CornerTableAdaptiveRefinement::makeUndo( )
{
    //Do nothing is the stack is empty.
    if (_undos.empty( ))
    {
        return;
    }

    //Get the top of the stack.
    UndoRedoRefinement* undo = _undos.top( );

    //Remove the top of the stack.
    _undos.pop( );

    //Perform the undo operations.
    while (undo->makeUndo( _triangleLevels, _vertexLevels ));

    //Add the redo operations on stack of redo.
    _redos.push( undo );
}



bool CornerTableAdaptiveRefinement::redoIsEmpty( )
{
    return _redos.empty( );
}



bool CornerTableAdaptiveRefinement::undoIsEmpty( )
{
    return _undos.empty( );
}
