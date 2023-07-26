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

#ifndef CORNERTABLEMULTILEVELSUBDIVISION_H
#define	CORNERTABLEMULTILEVELSUBDIVISION_H

#include "CornerTable.h"
#include "UndoRedoRefinement.h"
#include "ManagerTriangleLevels.h"
#include "ManagerUndoOperations.h"

#include <queue>
#include <vector>
#include <set>
#include <stack>


/**@class NodeOfSearch
 * Information to equeue during traversing the mesh to perform the refinement.
 */
class NodeOfSearch
{
public:
    /**
     * Corner of the next vertex to process in search during the refinement.
     */
    CornerType _vertex;

    /**
     * Opposite corner to the edge that needs to be flipped after the vertex be
     * processed.
     */
    CornerType _edge;
};


/**@class Triangle
 * Structure to store the information of triangles to be refined. This structure
 * is used on the preprocessing step and to sort the triangles by level.
 */
class Triangle
{
public:
    /**
     * Triangle index.
     */
    CornerType _index;

    /**
     * Current refinement level.
     */
    LevelType _level;

    /**
     * Overload of operator to be used on sort step.
     * @param t - triangle to be tested.
     * @return - true if the triangles are in correct order and false otherwise.
     */
    bool operator<( const Triangle& t ) const
    {
        if ( _level < t._level )
        {
            return true;
        }
        else if ( _level == t._level )
        {
            return _index < t._index;
        }
        return false;
    }

    /**
     * Overload of operator to determine if two triangles are equal or not.
     * @param t - triangle to be tested.
     * @return - true if the triangles are equal and false otherwise.
     */
    bool operator==( const Triangle& t ) const
    {
        return _index == t._index;
    }
};


class CornerTableAdaptiveRefinement
{
public:


    /**
     * Enumerator to represent the type of triangulation to be performed.
     */
    enum Triangulation
    {
        ARTMe = 0,
        REDGREEN = 1,
        INCREMENTAL = 2
    };

    /**
     * Constructor thar receives a surface as parameter.
     * @param surface - surface represented with corner table topological data
     * structure.
     */
    CornerTableAdaptiveRefinement( CornerTable *surface );

    /**
     * Read a mesh with level and sublevel of triangles and vertices.
     * @param path - path of the file.
     */
    CornerTableAdaptiveRefinement( const char* path );


    /**
     * Destructor.
     */
    virtual ~CornerTableAdaptiveRefinement( );

    /**
     * Get the refinement level of the triangle.
     * @param t - triangle index.
     * @return - refinement level of the triangle.
     */
    inline LevelType getTriangleLevel( const CornerType t ) const
    {
        return getLevel( &_triangleLevels[t] );
    }

    /**
     * Get refinement sublevel of the triangle.
     * @param t - triangle index.
     * @return - refinement sublevel of the triangle.
     */
    inline unsigned char getTriangleSubLevel( const CornerType t ) const
    {
        return getSubLevel( &_triangleLevels[t] );
    }

    /**
     * Gete the current surface.
     * @return - a pointer to corner table.
     */
    CornerTable *getSurface( );

    /**
     * Get the refinement level of the vertex.
     * @param v - vertex index.
     * @return - refinement level of the vertex.
     */
    LevelType getVertexLevel( const CornerType v );

    /**
     * Perform an undo step.
     */
    void makeUndo( );

    /**
     * Perform a redo step.
     */
    void makeRedo( );

    /**
     * Perform the global refinement on mesh. In this refinement, all triangles
     * are uniformly refined.
     * @param undo - true to enable the undo/redo mechanism and false otherwise.
     * @param computeNewVertex - callback to compute new vertex position.
     */
    void meshGlobalRefinement( const bool undo = false,
                               void (*computeNewVertex )( const CornerTable* cornerTable,
                               const CornerType corner, double *attributes ) = CornerTableAdaptiveRefinement::computeMiddlePoint );

    /**
     * Perform the local refinement on mesh. In this refinement, a subset of
     * triangles are refined. However, in a second step all triangles will be
     * considered as a coarse mesh, i.e., even that triangles that already began
     * to be refined will be refined as an initial triangle. 
     * @param trianglesForRefine - a set of triangles to be refined.
     * @param undo - true to enable the undo/redo mechanism and false otherwise.
     * @param computeNewVertex - callback to compute new vertex position.
     */
    void meshLocalRefinement( const std::set<CornerType>& trianglesForRefine,
                              const bool undo = false,
                              void (*computeNewVertex )( const CornerTable* cornerTable,
                              const CornerType corner, double *attributes ) = CornerTableAdaptiveRefinement::computeMiddlePoint );

    /**
     * Perform the adaptive refinement on a set of triangles. In this refinement,
     * a subset of triangles are refined and the algorithm is able to identify
     * triangles that already began to be refined and complete this refinement.
     * @param trianglesForRefine - a set of triangles to be refined.
     * @param undo - true to enable the undo/redo mechanism and false otherwise.
     * @param triangulation - the kind of triangulation to be generated in the
     * refined mesh. The possible types are: ARTMe, REDGREEN and INCREMENTAL.
     * @param ratio - the ratio used on incremental triangulation. This parameter
     * is used just in case of incremental triangulation.
     * @param computeNewVertex - callback to compute new vertex position.
     */
    void meshAdaptiveRefinement( const std::vector<CornerType>& trianglesForRefine,
                                 const bool undo = false, const Triangulation triangulation = ARTMe,
                                 const unsigned int ratio = 1,
                                 void (*computeNewVertex )( const CornerTable* cornerTable,
                                 const CornerType corner, double *attributes ) = CornerTableAdaptiveRefinement::computeMiddlePoint );

    /**
     * Perform the adaptive refinement on a set of triangles. In this refinement,
     * a subset of triangles are refined and the algorithm is able to identify
     * triangles that already began to be refined and complete this refinement.
     * @param trianglesForRefine - a set of triangles to be refined.
     * @param undo - true to enable the undo/redo mechanism and false otherwise.
     * @param triangulation - the kind of triangulation to be generated in the
     * refined mesh. The possible types are: ARTMe, REDGREEN and INCREMENTAL.
     * @param ratio - the ratio used on incremental triangulation. This parameter
     * is used just in case of incremental triangulation.
     * @param computeNewVertex - callback to compute new vertex position.
     */
    void meshAdaptiveRefinement( const std::set<CornerType>& trianglesForRefine,
                                 const bool undo = false, const Triangulation triangulation = ARTMe,
                                 const unsigned int ratio = 1,
                                 void (*computeNewVertex )( const CornerTable* cornerTable,
                                 const CornerType corner, double *attributes ) = CornerTableAdaptiveRefinement::computeMiddlePoint );

    /**
     * Return if the stack's redo is empty or not.
     * @return - true stack's redo is empty or false otherwise.
     */
    bool redoIsEmpty( );

    /**
     * Return if the stack's undo is empty or not.
     * @return - true stack's undo is empty or false otherwise.
     */
    bool undoIsEmpty( );

    /**
     * Write the current mesh on text file.
     * @param path - path of the file to be write.
     */
    void writeMesh( const char* path );

    /**
     * Writhe an OFF file.
     * @param path - path of the file to be write.
     */
    void writeMeshOFFFIle( const char* path );
private:

    /**
     * Pointer to surface stored using the corner table topological data
     * structure.
     */
    CornerTable *_surface;

    /**
     * Vetor utilizado apenas para evitar overhead de copias.
     */
    std::vector<CornerType> _cornerNeighbours;

    /**
     * Vector used to store the level and sublevel of the triangles.
     */
    std::vector<LevelType> _triangleLevels;

    /**
     * Vector used to store the level of the vertices.
     */
    std::vector<LevelType> _vertexLevels;

    /**
     * Stack to store the undo operations.
     */
    std::stack<UndoRedoRefinement*> _undos;

    /**
     * Stack to store the redo operations.
     */
    std::stack<UndoRedoRefinement*> _redos;

    /**
     * Auxiliary vector to be used by split operation to avoid relocations.
     */
    double *_coordinates;
private:
    /**
     * Perform the preprocessing step for the ARTMe refinement algorithm. This 
     * function verify is there is selected triangles to be refined with 
     * neighbors that have a lower refinement level. In this case, that neighbors
     * triangles are also selected to be refined in order to guarantee a balanced
     * mesh. Another verification is if triangles with sublevel 3 is selected to
     * be refined. In this case, the neighbor triangle with sublevel 2 or 1 
     * replace this triangle on the input.
     * @param trianglesForRefine - current selected triangles to be refined.
     */
    void preprocessingARTMe( std::set<Triangle>& trianglesForRefine );

    /**
     * Perform the preprocessing step in order to generate a red-green
     * triangulation. This function verify is there is selected triangles to be
     * refined with neighbors that have a lower refinement level. In this case,
     * that neighbors triangles are also selected to be refined in order to 
     * guarantee a balanced mesh. Besides that, this function select a triangle
     * with sublevel 1 or 2 to be refined always your mate will be refined in
     * order to avoid the generation of a triangle with sublevel 3.
     * @param trianglesForRefine - current selected triangles to be refined.
     */
    void preprocessingRedGreen( std::set<Triangle>& trianglesForRefine );

    /**
     * Perform the preprocessing step in order to generate an incremental
     * triangulation. This function expands the selected area add the functions
     * in a radius r to be refined.
     * @param trianglesForRefine - current selected triangles to be refined.
     * @param radius - radius r to be used on expand step.
     */
    void preprocessingIncremental( std::set<Triangle>& trianglesForRefine, unsigned int radius );

    /**
     * For each triangle with sublevel 3, replace it by the mate triangle with
     * sublevel 1 or 2.
     * @param triangle - index of the triangle with sublevel 3.
     * @param t - mate triangle to be filled.
     * @return - true is the mate triangle was found and false otherwise. The
     * expected return is always true.
     */
    bool addNeighbor3SublevelTriangles( CornerType triangle, Triangle &t );

    /**
     * Add neighbor triangles with lower level to be refined.
     * @param trianglesForRefine - current selected triangles to be refined.
     * @param newTriangles - vector that are being used to iterate on triangles
     * and must be updated on the function.
     * @param t - current triangle.
     */
    void addLowerLevelTriangles( std::set<Triangle>& trianglesForRefine,
                                 std::vector<Triangle>& newTriangles,
                                 const Triangle& t );

    /**
     * Add siblings triangles to be refined, i.e., whenever a triangle with
     * sublevel 1 or 2 will be refined, its sibling triangle with sublevel 2 or
     * 1, respectively, needs to be selected to be refined in order to generate
     * a red-green triangulation.
     * @param trianglesForRefine - current selected triangles to be refined.
     * @param newTriangles - vector that are being used to iterate on triangles
     * and must be updated on the function.
     * @param visitedTriangles - vector used to speed up the preprocessing step.
     * @param t - current triangle.
     */
    void addSiblingTriangles( std::set<Triangle>& trianglesForRefine,
                              std::vector<Triangle>& newTriangles,
                              std::vector<bool>& visitedTriangles,
                              const Triangle& t );

    /**
     * Add opposite trinagles that will generate triangles with sublevel 3, or
     * with two or more cracks (this triangles just can have sublevel 0) in 
     * order to generate a red-green triangulation.
     * @param trianglesForRefine - current selected triangles to be refined.
     * @param newTriangles - vector that are being used to iterate on triangles
     * and must be updated on the function.
     * @param visitedTriangles - vector used to speed up the preprocessing step.
     * @param t - current triangle.
     */
    void addOpositeTriangles( std::set<Triangle>& trianglesForRefine,
                              std::vector<Triangle>& newTriangles,
                              std::vector<bool>& visitedTriangles,
                              const Triangle& t );

    /**
     * From a set of vertices, compute the triangles formed by them in the
     * incremental preprocessing.
     * @param vertices - vertices to compute the triangles.
     * @param triangles - triangles formed by the vertices.
     */
    void computeTrianglesFormedByVertexs( const std::set<CornerType>& vertices,
                                          std::set<Triangle>& triangles );

    /**
     * Expand the selected area on preprocessing step in incremental preprocessing
     * in order to add that triangle that are in the r-ring neighborhood.
     * @param corners - initial corners of the selected area to begin the search.
     * @param radius - radius to search.
     * @param vertices - vertices on the r-ring neighborhood.
     */
    void getRRingFromVertices( const std::set<CornerType>& corners,
                               const unsigned int radius,
                               std::set<CornerType>& vertices );

    /**
     * Perform the refinement of the selected triangles after preprocessing step.
     * Perform the refinement of triangles from the lowest to the highest
     * refinement level avoiding T-vertices and any inconsistencies in the
     * refinement algorithm.
     * @param trianglesForRefine - set of triangles that must be refined.
     * @param undo - true to enable the undo/redo mechanism and false otherwise.
     * @param computeNewVertex - callback to compute new vertex position.
     * @return - return the stack of undo operations.
     */
    UndoRedoRefinement* makeRefinement( const std::set<Triangle>& trianglesForRefine,
                                        const bool undo,
                                        void (*computeNewVertex )( const CornerTable* cornerTable,
                                        const CornerType corner, double *attributes ) );

    /**
     * Perform the refinement of a single and connected component. The selected
     * triangles are implicit on processedFaces vector. Each selected triangle 
     * are marked as processed.
     * @param processedFaces - vector to mark if a triangles was processed or not.
     * @param initialCorner - corner to initialize the refinement process.
     * @param undo - true to enable the undo/redo mechanism and false otherwise.
     * @param computeNewVertex - callback to compute new vertex position.
     */
    void meshRefineTriangles( std::vector<bool>& processedTriangles,
                              const CornerType initialCorner, const bool undo,
                              UndoRedoRefinement* currentUndoRedo,
                              void (*computeNewVertex )( const CornerTable* cornerTable,
                              const CornerType corner, double *attributes ) );

    /**
     * Update the refinement level of the faces after an Edge Split operation.
     * This function just update the refinement level inside a single triangle
     * that is refined in two.
     * @param t1 - index of the original triangle.
     * @param t2 - index of the new created triangle.
     */
    void processTriangleLevel( CornerType t1, CornerType t2 );

    /**
     * Verify if the triangles was processed and mark them as processed if it is
     * needs. Triangles with sublevel 0 or 3 must be marked as processed, once
     * this update must be performed together with the split operation.
     * @param processedTriangles - vector to mark if a triangles was processed
     * or not.
     * @param t1 - index of the original triangle.
     * @param t2 - index of the new created triangle.
     */
    void processTriangles( std::vector<bool>& processedTriangles, CornerType t1, CornerType t2 );

    /**
     * Perform the Edge Split operation and update the triangle status in the
     * local and global refinement.
     * @param processedTriangles - vector to mark if a triangles was processed
     * or not. This vector is updated inside the function.
     * @param visitedTriangles - vector to mark if a triangles was visited
     * or not. This vector is updated inside the function.
     * @param cornersForSplit - queue of elements to refinement process. This
     * queue has a pair of vertex to be processed and edge to be flipped.
     * @param cornerOperation - index corner to perform the edge split operation.
     * @param undoRedo - object to store the operations to perform undo.
     * @param undo - true to enable the undo/redo mechanism and false otherwise.
     * @param computeNewVertex - callback to compute new vertex position.
     */
    void runSplit( std::vector<bool>& processedTriangles,
                   std::vector<bool>& visitedTriangles,
                   std::queue<NodeOfSearch>& cornersForSplit,
                   const CornerType cornerOperation,
                   UndoRedoRefinement* undoRedo, const bool undo,
                   void (*computeNewVertex )( const CornerTable* cornerTable,
                   const CornerType corner, double *attributes ) );

    /**
     * Refine a single and connected component triangles in the global and local
     * refinement.
     * @param t - index of initial triangle.
     * @param processedTriangles - vector to mark if a triangles was processed
     * or not. This vector is updated inside the function.
     * @param visitedTriangles - vector to mark if a triangles was visited
     * or not. This vector is updated inside the function.
     * @param cornersForSplit - queue of elements to refinement process. This
     * queue has a pair of vertex to be processed and edge to be flipped.
     * @param undoRedo - object to store the operations to perform undo.
     * @param undo - true to enable the undo/redo mechanism and false otherwise.
     * @param computeNewVertex - callback to compute new vertex position.
     */
    void refineTriangles( const CornerType t,
                          std::vector<bool>& processedTriangles,
                          std::vector<bool>& visitedTriangles,
                          std::queue<NodeOfSearch>& cornersForSplit,
                          UndoRedoRefinement* undoRedo, const bool undo,
                          void (*computeNewVertex )( const CornerTable* cornerTable,
                          const CornerType corner, double *attributes ) );

    /**
     * Perform the split operation to adaptive refinement and update the triangles
     * status.
     * @param processedFaces - vector to mark if a triangles was processed
     * or not. This vector is updated inside the function.
     * @param cornersForSplit - queue of elements to refinement process. This
     * queue has a pair of vertex to be processed and edge to be flipped.
     * @param cornerOperation - index corner to perform the edge split operation.
     * @param undoRedo - object to store the operations to perform undo.
     * @param undo - true to enable the undo/redo mechanism and false otherwise.
     * @param computeNewVertex - callback to compute new vertex position.
     */
    void runSplitAdaptiveRefinement( std::vector<bool>& processedTriangles,
                                     std::queue<NodeOfSearch>& cornersForSplit,
                                     const CornerType cornerOperation,
                                     UndoRedoRefinement* undoRedo, const bool undo,
                                     void (*computeNewVertex )( const CornerTable* cornerTable,
                                     const CornerType corner, double *attributes ) );

    /**
     * Try to flip an edge on the local and global refinement.
     * @param corner - opposite corner to the edge.
     * @param undo - true to enable the undo/redo mechanism and false otherwise.
     * @param undoRedo - object to store the operations to perform undo.
     */
    void tryFlipEdge( const CornerType corner, bool undo, UndoRedoRefinement* undoRedo );

    /**
     * Try to flip an edge if it is possible, i.e., the edge has two triangles
     * with sublevel 3.
     * @param corner - opposite corner to the edge.
     * @param undo - true to enable the undo/redo mechanism and false otherwise.
     * @param undoRedo - object to store the operations to perform undo.
     */
    void tryFlipEdgeAdaptiveRefinement( const CornerType corner, bool undo,
                                        UndoRedoRefinement* undoRedo );

    /**
     * Update the stack of undo and delete the stack of redo.
     * @param undo - undo stack.
     */
    void updateUndos( const UndoRedoRefinement* undo );

    /**
     * Function to compute the point at the opposite edge. This points is calculated
     * in the middle of the edge.
     * @param cornerTable - surface represented as a corner table.
     * @param c - index of corner opposite to the edge.
     * @param attributes - vector that must be filled with the properties of the
     * new vertex.
     */
    static void computeMiddlePoint( const CornerTable* cornerTable,
                                    const CornerType c, double* attributes );
};

#endif	/* CORNERTABLEMULTILEVELSUBDIVISION_H */

