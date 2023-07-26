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

#ifndef DEFINITIOTYPES_H
#define	DEFINITIOTYPES_H

/**
 * Type to represent the corners index, triangle index and vertex index. If the
 * undo/redo mechanism is used, the representation lose two bits.
 */
typedef unsigned int CornerType;

/**
 * Type to represent the level and sublevel of triangles and the level of the
 * vertices.
 */
typedef unsigned char LevelType;


#endif	/* DEFINITIOTYPES_H */

