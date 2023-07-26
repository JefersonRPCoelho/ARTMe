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
#include "ManagerTriangleLevels.h"



void setTriangleLevels( LevelType *triangle, LevelType level, LevelType sublevel )
{
    (*triangle) = (level << 2) | sublevel;
}



void setTriangleSubLevel( LevelType *triangle, LevelType sublevel )
{
    (*triangle) = (((*triangle) & ~0x3) | (sublevel));
}



LevelType getLevel( const LevelType *triangle )
{
    return (*triangle) >> 2;
}



LevelType getSubLevel( const LevelType *triangle )
{
    return (*triangle)& 0x3;
}



void decrementSubLevel( LevelType *triangle )
{
    LevelType sublevel = getSubLevel( triangle );

    LevelType level = getLevel( triangle );

    level += (sublevel + 3) / 4 - 1;
    sublevel = (sublevel + 3) % 4;

    setTriangleLevels( triangle, level, sublevel );
}



void incrementSubLevel( LevelType *triangle, LevelType increment )
{
    LevelType sublevel = getSubLevel( triangle );

    LevelType level = getLevel( triangle );

    level += (sublevel + increment) / 4;
    sublevel = (sublevel + increment) % 4;

    setTriangleLevels( triangle, level, sublevel );
}



