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


#ifndef UNDO_REDO_H
#define	UNDO_REDO_H

#include <stack>
#include "CornerTable.h"


/**
 * Class to store the undo/redo information and execute the operations.
 */

class UndoRedoRefinement
{
public:


    /**
     * Possible operations.
     */
    enum OperationType
    {
        WELD = 0,
        FLIP = 1,
        SPLIT = 2,
        UNFLIP = 3
    };

    /**
     * Constructor that receives a surface.
     */
    UndoRedoRefinement( CornerTable* surface );

    /**
     * Copy construtor.
     * @param orig - an object.
     */
    UndoRedoRefinement( const UndoRedoRefinement& orig );

    /**
     * Perform the redo.
     * @param triangleLevels - vector with the triangle refinement level to be
     * updated inside the function.
     * @param vertexLevels - vector with the vertex refinement level to be
     * updated inside the function.
     * @return - true is it is possible to execute the redo.
     */
    bool makeRedo( std::vector<LevelType>& triangleLevels,
                   std::vector<LevelType>& vertexLevels );

    /**
     * Perform the undo.
     * @param triangleLevels - vector with the triangle refinement level to be
     * updated inside the function.
     * @param vertexLevels - vector with the vertex refinement level. In this case
     * the size of vector is decremented.
     * @return - true is it is possible to execute the undo.
     */
    bool makeUndo( std::vector<LevelType>& triangleLevels,
                   std::vector<LevelType>& vertexLevels );

    /**
     * Define new information to execute undo operation.
     * @param operation - an integer with the operation and the corner to apply
     * the undo operation.
     */
    void addUndo( const CornerType operation );

    /**
     * Clear the undo and redo stack.
     */
    void reset( );

    /**
     * Destructor.
     */
    virtual ~UndoRedoRefinement( );

    /*
     * Overload of assignment operator.
     */
    const UndoRedoRefinement& operator=( const UndoRedoRefinement& orig );

    /**
     *Verify if the both stacks are empty.
     * @return - true if the two stacks are empty and false otherwise.
     */
    bool isEmpty( ) const;

    /**
     * Update the refinement level of the faces after an Edge Split operation.
     * This function just update the refinement level inside a single triangle
     * that is refined in two.
     * @param triangleLevels - vector with the triangle refinement level to be
     * updated inside the function.
     * @param t1 - index of the original triangle.
     * @param t2 - index of the new created triangle.
     */
    void processTriangleLevel( std::vector<LevelType>& triangleLevels, CornerType f1, CornerType f2 );
private:

    /**
     * Surface represented with the corner table topological data structure.
     */
    CornerTable* _surface;

    /**
     * Undo stack.
     */
    std::stack<CornerType> _undo;

    /**
     * Redo stack.
     */
    std::stack<CornerType> _redo;
};

#endif	/* UNDO_REDO_H */
