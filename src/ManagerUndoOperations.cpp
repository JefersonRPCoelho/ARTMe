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

#include <stdlib.h>
#include <stdio.h>
#include "ManagerUndoOperations.h"



void setOperation( CornerType *operation, CornerType corner, CornerType operationType )
{
    ( *operation ) = ( corner << 2 ) | operationType;
}



void setOperationType( CornerType *operation, CornerType operationType )
{
    ( *operation ) = ( ( *operation ) & ~0x3 ) | ( operationType );
}



void setCornerOperation( CornerType *operation, CornerType corner )
{
    ( *operation ) = ( corner << 2 ) | ( ( *operation )& 0x3 );
}



CornerType getCornerOperation( const CornerType *operation )
{
    return (*operation ) >> 2;
}



unsigned char getOperationType( const CornerType *operation )
{
    return (*operation )& 0x3;
}


