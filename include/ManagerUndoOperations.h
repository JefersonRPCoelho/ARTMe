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


#ifndef MANAGERUNDOOPERATIONS_H
#define	MANAGERUNDOOPERATIONS_H
#include  "DefinitionTypes.h"
#ifdef __cplusplus
extern "C"
{
#endif


    /**
     * Define the corner and the operation type.
     * @param operation - variable to store the informations of undo/redo
     * mechanism.
     * @param corner - corner where the operation must be executed.
     * @param operationType - operation type that must be executed.
     */
    void setOperation( CornerType *operation, CornerType corner, CornerType operationType );

    /**
     * Define the operation type.
     * @param operation - variable to store the informations of undo/redo
     * mechanism.
     * @param operationType - operation type that must be executed.
     */
    void setOperationType( CornerType *operation, CornerType operationType );

    /**
     * Define the corner to apply the operation.
     * @param operation - variable to store the informations of undo/redo
     * mechanism.
     * @param corner - corner where the operation must be executed.
     */
    void setCornerOperation( CornerType *operation, CornerType corner );

    /**
     * Get the corner to apply the operation.
     * @param operation - variable to store the informations of undo/redo
     * mechanism.
     * @return - corner to apply the operation.
     */
    CornerType getCornerOperation( const CornerType *operation );

    /**
     * Get the operation type.
     * @param operation operation - variable to store the informations of undo/redo
     * mechanism.
     * @return - operation type.
     */
    unsigned char getOperationType( const CornerType *operation );

#ifdef __cplusplus
}
#endif


#endif	/* MANAGERUNDOOPERATIONS_H */

