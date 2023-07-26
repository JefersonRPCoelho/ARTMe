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


#ifndef MANAGERFACELEVELS_H
#define	MANAGERFACELEVELS_H
#include "DefinitionTypes.h"
#ifdef __cplusplus
extern "C"
{
#endif


    /**
     * Define the level and sublevel of a triangle.
     * @param triangle - variable that represent the level of sublevel of a
     * triangle.
     * @param level - new level of the triangle.
     * @param sublevel - new sublevel of the triangle.
     */
    void setTriangleLevels( LevelType *triangle, LevelType level, LevelType sublevel );

    /**
     * Define the sublevel of the triangle.
     * @param triangle - variable that represent the level of sublevel of a
     * triangle.
     * @param sublevel - new sublevel of the triangle.
     */
    void setTriangleSubLevel( LevelType *triangle, LevelType sublevel );

    /**
     * Get the triangle refinement level.
     * @param triangle - variable that represent the level of sublevel of a
     * triangle.
     * @return - triangle refinement level.
     */
    LevelType getLevel( const LevelType *triangle );

    /**
     * Get the triangle refinement sublevel.
     * @param triangle - variable that represent the level of sublevel of a
     * triangle.
     * @return - triangle refinement sublevel.
     */
    LevelType getSubLevel( const LevelType *triangle );

    /**
     * Decrement the sublevel of the triangle.
     * @param triangle - variable that represent the level of sublevel of a
     * triangle.
     */
    void decrementSubLevel( LevelType *triangle );

    /**
     * Increment the sublevel of the triangle.
     * @param triangle - variable that represent the level of sublevel of a
     * triangle.
     * @param increment - increment to sum in the triangle refinement sublevel.
     */
    void incrementSubLevel( LevelType *triangle, LevelType increment = 1 );

#ifdef __cplusplus
}
#endif


#endif	/* MANAGERFACELEVELS_H */

