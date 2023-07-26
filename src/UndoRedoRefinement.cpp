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


#include "UndoRedoRefinement.h"
#include "ManagerTriangleLevels.h"
#include "ManagerUndoOperations.h"
#include <cstdio>



UndoRedoRefinement::UndoRedoRefinement( CornerTable* surface )
{
    _surface = surface;
}



UndoRedoRefinement::~UndoRedoRefinement( )
{
}



UndoRedoRefinement::UndoRedoRefinement( const UndoRedoRefinement& orig )
{
    this->_redo = orig._redo;
    this->_undo = orig._undo;
    this->_surface = orig._surface;
}



void UndoRedoRefinement::processTriangleLevel( std::vector<LevelType>& triangleLevels, CornerType t1, CornerType t2 )
{
    //Get the sublevel of the original triangle.
    const LevelType sublevel = getSubLevel( &triangleLevels[t1] );

    //Set the refinement level information of the new triangle equal to the old
    //triangle.
    triangleLevels.push_back( triangleLevels[t1] );

    switch (sublevel)
    {
        case 0:
        {
            //For sublevel 0, the old face will have sublevel 1 and the new face
            //sublevel 2.
            incrementSubLevel( &triangleLevels[t1] );
            incrementSubLevel( &triangleLevels[t2], 2 );
        }
            break;
        case 1:
        {
            //For sublevel 1, the old face will have the sublevel 3 and the new
            //face will have sublevel 0 in the next refinement level.
            incrementSubLevel( &triangleLevels[t1], 2 );
            incrementSubLevel( &triangleLevels[t2], 3 );
        }
            break;
        case 2:
        {
            //For sublevel 2, the oold face will have sublevel 0 in the next
            //refinement level and the new face will have sublevel 3.
            incrementSubLevel( &triangleLevels[t1], 2 );
            incrementSubLevel( &triangleLevels[t2] );
        }
            break;
        default:
        {
            printf( "Error: trying to process an invalid sublevel in Edge Split Operation\n" );
        }
    }
}



bool UndoRedoRefinement::makeRedo( std::vector<LevelType>& faceLevels,
                                   std::vector<LevelType>& vertexLevels )
{
    //Do nothing if the redo stack is empty.
    if (_redo.size( ) == 0)
    {
        return false;
    }

    //Get the operation on the top of the stack.
    CornerType op = _redo.top( );

    //Remove the top element of the stack.
    _redo.pop( );

    //Create a variable to represent the undo operation..
    CornerType undo;

    switch (getOperationType( &op ))
    {
        case FLIP:
        {
            CornerType corner = getCornerOperation( &op );

            //Get the involved triangles.
            CornerType oppositeTriangle = _surface->cornerOpposite( corner ) / 3;

            //Decrement the subnivel of the triangles..
            incrementSubLevel( &faceLevels[oppositeTriangle] );
            incrementSubLevel( &faceLevels[corner / 3] );

            //Perform the Edge Flip operation.
            _surface->edgeFlip( corner );

            //Create an undo to flip the edge
            CornerType newCorner = _surface->cornerNext( corner );
            setOperation( &undo, newCorner, UNFLIP );
            _undo.push( undo );
            break;
        }
        case SPLIT:
        {
            //Get the number of coordinates by vertex.
            const unsigned int n = _surface->getNumberAttributesByVertex( );

            //Allocate memory to store the new vertex.
            double * attributes = new double[n];

            //Get the geometry information of the  surface.
            const double* attributesCornerTable = _surface->getAttributes( );

            //Get the index of the new vertex.
            const CornerType vertex = _surface->getNumberVertices( );

            //Copy the vertex information.
            for (unsigned int i = 0; i < n; i++)
            {
                attributes[i] = attributesCornerTable[n * vertex + i];
            }

            CornerType corner = getCornerOperation( &op );

            //Get the opposite corner.
            CornerType opposite = _surface->cornerOpposite( corner );

            //Get the number of triangles on mesh.
            const CornerType numTriangles = _surface->getNumTriangles( );

            //Set the refinement level of the new vertex.
            vertexLevels.push_back( getLevel( &faceLevels[corner / 3] ) + 1 );

            //Process the triangle level on current triangle.
            processTriangleLevel( faceLevels, corner / 3, numTriangles );

            //Process the triangle level on opposite triangle.
            if (opposite != CornerTable::BORDER_CORNER)
            {
                processTriangleLevel( faceLevels, opposite / 3, numTriangles + 1 );
            }

            //Perform the split operation.
            _surface->edgeSplit( corner, attributes );

            //Delete auxiliary vector.
            delete []attributes;

            //Create an undo to weld the edge
            CornerType newCorner = _surface->cornerNext( corner );
            setOperation( &undo, newCorner, WELD );
            _undo.push( undo );
            break;
        }
        default:
            printf( "Error performing undo operation\n" );
            return false;
            break;
    }
    return true;
}



bool UndoRedoRefinement::makeUndo( std::vector<LevelType>& triangleLevels,
                                   std::vector<LevelType>& vertexLevels )
{
    //Do nothing if the redo stack is empty.
    if (_undo.size( ) == 0)
    {
        return false;
    }

    //Get the operation on the top of the stack.
    CornerType op = _undo.top( );

    //Remove the top element of the stack.
    _undo.pop( );

    //Create a variable to represent the undo operation..
    CornerType redo;

    switch (getOperationType( &op ))
    {
        case WELD:
        {
            CornerType corner = getCornerOperation( &op );

            //Triangle refinement sublevel after weld operation.
            LevelType triangleSubLevel = 0, oppositeTriangleSubLevel = 0;

            //Index of the current triangle..
            CornerType triangle = corner / 3;

            //Get the left triangle..
            CornerType oppositeTriangle = _surface->cornerLeft( corner );
            if (oppositeTriangle != CornerTable::BORDER_CORNER)
            {
                oppositeTriangle /= 3;
            }

            //Get the sublevel refinement of current triangle..
            LevelType triangleWeldSubLevel = getSubLevel( &triangleLevels[triangle] );

            //Get the sublevel after weld operation.
            if (triangleWeldSubLevel == 1)
            {
                triangleSubLevel = 0;
            }
            else if (triangleWeldSubLevel == 3)
            {
                triangleSubLevel = 1;
            }
            else if (triangleWeldSubLevel == 0)
            {
                triangleSubLevel = 2;
            }
            else
            {
                printf( "Erro no weld\n" );
            }

            //Get the sublevel of the opposite triangle.
            if (oppositeTriangle != CornerTable::BORDER_CORNER)
            {
                LevelType oppositeFaceWeldSubLevel = getSubLevel( &triangleLevels[oppositeTriangle] );
                if (oppositeFaceWeldSubLevel == 2)
                {
                    oppositeTriangleSubLevel = 0;
                }
                else if (oppositeFaceWeldSubLevel == 0)
                {
                    oppositeTriangleSubLevel = 1;
                }
                else if (oppositeFaceWeldSubLevel == 3)
                {
                    oppositeTriangleSubLevel = 2;
                }
                else
                {
                    printf( "Error processing weld on undo operation\n" );
                }
            }

            CornerType newCorner = _surface->cornerPrevious( corner );
            setOperation( &redo, newCorner, SPLIT );
            _redo.push( redo );

            //Perform weld operation..
            _surface->edgeWeld( corner );

            //Define the corrent level and sublevel refinement,
            triangle = newCorner / 3;
            oppositeTriangle = _surface->cornerOpposite( newCorner );

            //Decrement the sublevel of the current triangle.
            decrementSubLevel( &triangleLevels[triangle] );

            //Define the sublevel of the triangle.
            setTriangleSubLevel( &triangleLevels[triangle], triangleSubLevel );

            if (oppositeTriangle != CornerTable::BORDER_CORNER)
            {
                oppositeTriangle /= 3;

                //Decrement the sublevel of the current triangle.
                decrementSubLevel( &triangleLevels[oppositeTriangle ] );

                //Define the sublevel of the triangle.
                setTriangleSubLevel( &triangleLevels[oppositeTriangle ], oppositeTriangleSubLevel );

                //Remove a triangle in the final of the vector.
                triangleLevels.resize( triangleLevels.size( ) - 1 );
            }

            //Remove a triangle in the final of the vector.
            triangleLevels.resize( triangleLevels.size( ) - 1 );

            //Remove a vertex in the final of the vector.
            vertexLevels.resize( vertexLevels.size( ) - 1 );

            break;
        }

        case UNFLIP:
        {
            CornerType corner = getCornerOperation( &op );

            //Get the triangles involved in the flip operation.
            CornerType oppositeTriangle = _surface->cornerOpposite( corner ) / 3;

            //Decrement the sublevel of the both triangles.
            decrementSubLevel( &triangleLevels[oppositeTriangle] );
            decrementSubLevel( &triangleLevels[corner / 3] );

            //Perform the unflip operation.
            _surface->edgeUnflip( corner );

            CornerType newCorner = _surface->cornerPrevious( corner );
            setOperation( &redo, newCorner, FLIP );
            _redo.push( redo );
            break;
        }
        default:
            printf( "Error performing undo operation\n" );
            return false;
            break;
    }
    return true;
}



void UndoRedoRefinement::reset( )
{
    while (!_redo.empty( ))
    {
        _redo.pop( );
    }
    while (!_undo.empty( ))
    {
        _undo.pop( );
    }
}



void UndoRedoRefinement::addUndo( const CornerType operation )
{
    _undo.push( operation );
}



const UndoRedoRefinement& UndoRedoRefinement::operator=( const UndoRedoRefinement& orig )
{
    this->_redo = orig._redo;
    this->_undo = orig._undo;
    this->_surface = orig._surface;
    return *this;
}



bool UndoRedoRefinement::isEmpty( ) const
{
    return (_undo.empty( ) && _redo.empty( ) );
}
