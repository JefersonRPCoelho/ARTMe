# ARTMe - Adaptive Refinement for Triangle Mesh

## Authors
* Jéferson Coêlho
* [Marcelo Gattass](http://webserver2.tecgraf.puc-rio.br/~mgattass/)
* [Hélio Lopes](http://www-di.inf.puc-rio.br/~lopes//)

## ARTMe description

#### Summary
[ARTMe](https://link.springer.com/article/10.1007/s00366-018-0579-5) is an efficient array-based algorithm for adaptive triangle mesh refinement capable of interactively generating millions of triangles. The refinement algorithm satisfies important topological mesh properties, e.g., vertex valence control and a good mesh gradation. Furthermore, all local topological modifications of the triangle mesh are based on stellar operators implemented on top of the Corner-Table topological data structure, a well known memory-efficient data structure. This implementation provides a good balance in the trade-off between memory and processing time with support to undo and redo operations.

The refinement algorithm can also continue the triangle refinement process from the last performed operation. More precisely, once a triangle has been partially refined in a previous step, the algorithm can identify and perform the remaining sequence of steps necessary to complete the triangular refinement. This property classifies the proposed refinement algorithm as an adaptive one. This characteristic is especially important for iterative methods, with which we can refine the mesh step by step. All local topological modifications to the triangle mesh are based on the [stellar operators](https://doi.org/10.1016/j.cagd.2009.08.004).

The algorithm is also flexible. With few or no modifications, it can refine meshes similar to traditional algorithms, e.g., [simple triangulation](https://link.springer.com/chapter/10.1007/978-3-642-55787-3_19), Red-Green triangulation and [incremental triangulation](http://www.inderscience.com/offer.php?id=14467). It uses six references per triangle to represent the mesh connectivity, which is 2.75 times less than that of the Half-Edge. Although the Corner-Table uses vectors to store the topological information, the implementation was projected so as to eliminate the necessity of garbage collection.

#### Strength
ARTMe has the following good properties
* ARTMe is a **very efficient** adaptive refinement algorithm. In the worst case, the ARTMe algorithm generates more than 1 million triangles by second, even on a simple computer (Intel(R) Core(TM)2 Quad processor with 2.66-GHz CPUs
and 4 GB of 800-MHz DDR2 memory).
* ARTMe is **versatile** enough to refine triangle meshes with the same pattern of standard algorithms, like Simple Triangulation, Red Green Triangulation and Incremental triangulation.
* Although the literature reports that the "Red Green triangulation method is neither efficient nor simple", the ARTMe algorithm is capable to perform the Red Green triangulation very efficiently.
* Once ARTMe uses only *stellar operators*, it is  mathematically guaranteed that the **surface topology does not change** during the refinement process.
* ARTMe has the desirable property of generating **new vertices with valence six**, while the valences of the **old vertices do not change**. On the border of the refined area, the valence of a vertex can increase at most double, before returning to the initial value.
* Using the Corner-Table to replace the standard Half-Edge data structure, ARTMe **reduce the memory consumption** of the required topological information for mesh refinement by a factor of 2.75.
* ARTMe produces a mesh with **good mesh gradation**. In the resulting mesh the difference in refinement levels of two adjacent faces is not greater than one.
* At the same conditions, ARTMe triangulation pattern **generates fewer triangles** to satisfy the constraints, like a maximum threshold to an approximation.
* ARTMe **supports the undo and redo mechanism** without needing garbage collection operations.

#### Weakness
ARTMe has the following weakness
* ARTMe is not parallelizable.
* ARTMe increases the number of triangles by a factor of four.

You can find more information about ARTMe method [here](https://link.springer.com/article/10.1007/s00366-018-0579-5).

#### Citation
If you find this work helpful in your research, please cite:
```
@Article{Coelho2018,
    author   = "Co{\^e}lho, J{\'e}ferson and Gattass, Marcelo and Lopes, H{\'e}lio",
    title    = "ARTMe: a new array-based algorithm for Adaptive Refinement of Triangle Meshes",
    journal  = "Engineering with Computers",
    year     = "2018",
    month    = "Jan",
    day      = "18",
    abstract = "This work presents a new efficient array-based algorithm for adaptive mesh refinement capable of interactively generating millions of triangles. The new refinement algorithm satisfies important topological mesh properties, e.g., vertex valence control and a good mesh gradation. Furthermore, all local topological modifications of the triangle mesh are based on Stellar operators implemented on top of the Corner-Table topological data structure. This paper also shows that the proposed implementation provides a good balance in the trade-off between memory and processing time.",
    issn     = "1435-5663",
    doi      = "10.1007/s00366-018-0579-5",
    url      = "https://doi.org/10.1007/s00366-018-0579-5"
}
```

#### License
This software is provided by the [Tecgraf/PUC-Rio Institute](https://www.tecgraf.puc-rio.br/) "as is" and any express or implied warranties, including, but not limited to, the implied warranties of merchantability and fitness for a particular purpose are disclaimed. In no event shall the scalable software infrastructure project be liable for any direct, indirect, incidental, special, exemplary, or consequential damages (including, but not limited to, procurement of substitute goods or services; loss of use, data, or profits; or business interruption) however caused and on any theory of liability, whether in contract, strict liability, or tort (including negligence or otherwise) arising in any way out of the use of this software, even if advised of the possibility of such damage.

## Getting Started


These instructions will get you a copy of the project up and running on your local machine for development and testing purposes.

#### Prerequisites

This project has no dependencies.

#### Getting the project
To get the project, use
```
git clone  https://github.com/JefersonRPCoelho/ARTMe.git
```

#### Compiling

To compile this project and generate shared and static libraries you can use the available CMakeList:

From the project directory, create a build directory

```
mkdir build
cd build
```

Execute cmake

```
cmake ..
```
Compile the source code
```
make
```
Install the resulting binaries

```
make install
```

At this point, there is an executable file called *ARTMe* on *bin* directory and the library files *libARTMe.so* and *libARTMe.a* are on *lib* directory.
## Running

#### Running the tests

To run the resulting program, enter into *bin* directory and use

```
./ARTMe
```

#### Results

The sample program generates two refined models as example. The first one is a regular hexagon mesh that is refined using the ARTMe, Red-Green and Incremental triangulations. The second one refines a regular icosahedron in order to generate a sphere. The refinement process is performed on Hexagon' and Sphere' bottom and top using all three triangulation patterns.

The sample program generates files named as the follow pattern
```
<example><Triangulation><Version>.off
```
with
```
example = {hexagon, sphere}
Triangulation = {ARTMe, Incremental, RedGreen}
Version = {original, refined, originalCopy, refinedCopy}
```
The input to two exemples is an *original* mesh. After the execution of ARTMe algorithm a *refined* mesh is created. The undo mechanism is used over the *refined* mesh to generate a *original copy* mesh. Finally, the redo mechanism is used to generate the *refined copy* mesh.

Note that, the *original copy* and *refined copy* are completely identical to *original* and *refined* meshes, respectively, even in the index triangulation level and vertex coordinates.

##### Hexagon resulting meshes
| **ARTMe**  |  **Red-Green** | **Incremental**|
| :-------------: | :-------------: | :-------------: |
| ![Alt ARTMe](/images/hexARTMe.png)| ![Alt Red-Green](/images/hexRed-Green.png)| ![Alt incremental](/images/hexIncremental.png)|
| 1,388 triangles  |  1,476 triangles | 2,724 triangles|

##### Sphere resulting meshes
|                  **ARTMe**               |                       **Red-Green**              |                       **Incremental**                |
|                 :---:                |                 :-------------:              |                      :-------------:             |
| ![Alt ARTMe](/images/sphereARTMe.png)| ![Alt Red-Green](/images/sphereRed-Green.png)| ![Alt Incremental](/images/sphereIncremental.png)|
|                7,496 triangles       |                 7,632 triangles              |                      12,232 triangles            |




## ARTMe Library Interface
By including the file
```
CornerTableAdaptiveRefinement.h
```
you can construct an object to perform ARTMe algorithm. To do this, you need to store the triangle mesh in a *CornerTable* data structure and use the follow constructor:
```cpp
/**
  * Constructor that receives a surface as parameter.
  * @param surface - surface represented with corner table topological data
  * structure.
  */
CornerTableAdaptiveRefinement( CornerTable *surface );
```

#### Adaptive Refinement
Once you have an object, the library provides two functions to execute an adaptive refinement:
```cpp
/**
  * Perform the adaptive refinement on a set of triangles. In this refinement,
  * a subset of triangles are refined and the algorithm is able to identify
  * triangles that already began to be refined and complete this refinement.
  * @param trianglesForRefine - a set of triangles to be refined.
  * @param undo - true to enable the undo/redo mechanism and false otherwise.
  * @param triangulation - the kind of triangulation to be generated in the
  * refined mesh. The possible types are: ARTMe, REDGREEN and INCREMENTAL.
  * @param radius - the radius used on incremental triangulation. This parameter
  * is used just in case of incremental triangulation.
  * @param computeNewVertex - callback to compute new vertex position.
  */
void meshAdaptiveRefinement( const std::vector<CornerType>& trianglesForRefine,
                             const bool undo = false,
                             const Triangulation triangulation = ARTMe,
                             const unsigned int radius = 1,
                             void (*computeNewVertex )( const CornerTable *cornerTable,
                             const CornerType corner, double *attributes ) =
                             CornerTableAdaptiveRefinement::computeMiddlePoint );

/**
  * Perform the adaptive refinement on a set of triangles. In this refinement,
  * a subset of triangles are refined and the algorithm is able to identify
  * triangles that already began to be refined and complete this refinement.
  * @param trianglesForRefine - a set of triangles to be refined.
  * @param undo - true to enable the undo/redo mechanism and false otherwise.
  * @param triangulation - the kind of triangulation to be generated in the
  * refined mesh. The possible types are: ARTMe, REDGREEN and INCREMENTAL.
  * @param radius - the radius used on incremental triangulation. This parameter
  * is used just in case of incremental triangulation.
  * @param computeNewVertex - callback to compute new vertex position.
  */
void meshAdaptiveRefinement( const std::set<CornerType>& trianglesForRefine,
                             const bool undo = false,
                             const Triangulation triangulation = ARTMe,
                             const unsigned int radius = 1,
                             void (*computeNewVertex )( const CornerTable *cornerTable,
                             const CornerType corner, double *attributes ) =
                             CornerTableAdaptiveRefinement::computeMiddlePoint );
```

Note that the only difference between the functions is that the first one receives the list of triangles, that must to be refined, as a std::vector and the second one receives the same list as a std::set.

The only parameter required is the triangles list that will be refined, the remaining ones are optional.

By using all parameters, you can control
* *undo* - if the algorithm will store information about undo/redo operations. In case you do not explicitly set this parameter as **true**, the undo/redo mechanism will not be available to be used.
* *triangulation* - the triangulation pattern. ARTMe, REDGEEN or INCREMENTAL are the possible patterns.
* *radius* - the radius used by Incremental Triangulation. This parameter is ignored if *triangulation* != INCREMENTAL.
* *computeNewVertex* - a function to compute the new vertex coordinates.

The code
```cpp
....
CornerTableAdaptiveRefinement *r = CornerTableAdaptiveRefinement( s );
...
s->meshAdaptiveRefinement( t );
```
will refine all triangles on *t* list using ARTMe triangulation pattern. The refinement process will not store information about undo/redo operations and new vertices will be positioned at the midpoint of the edges.

To generate a sphere from an icosahedron, as the above examples, you can call the refinement function as follows
```cpp
s->meshAdaptiveRefinement( t, true, CornerTableAdaptiveRefinement::ARTMe,
                           0, computeNewVertexOverSphere );
```
The function *computeNewVertexOverSphere* is callback implemented to project each new point on sphere surface. As the icosahedron is centered at origin, it's enough to normalize the vector that represents the new vertex
```cpp
/**
  * Callback to compute the new vertices coordinates. In this case, the
  * callback set new coordinates in the unit sphere surface.
  * @param surface - current surface represented in the corner table data
  * structure.
  * @param corner - corner opposite to the edge where the new vertex will
  * be created.
  * @param attributes - vector to be set with the new coordinates.
  */
static void computeNewVertexOverSphere( const CornerTable *surface, const CornerType corner,
                                        double *attributes )
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

    //As the sphere has its center in (0, 0, 0), the solution is just compute the unit
    //vector to the middle point.
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
```

#### Global Refinement
If you intend to refine the entire mesh, it is more efficient to use the following function
```cpp
/**
  * Perform the global refinement on mesh. In this refinement, all triangles
  * are uniformly refined.
  * @param undo - true to enable the undo/redo mechanism and false otherwise.
  * @param computeNewVertex - callback to compute new vertex position.
  */
void meshGlobalRefinement( const bool undo = false,
                           void (*computeNewVertex )( const CornerTable *cornerTable,
                           const CornerType corner, double *attributes ) =
                           CornerTableAdaptiveRefinement::computeMiddlePoint );
```
The *meshGlobalRefinement* function implements an special case of the ARTMe algorithm. This function is faster than *meshAdaptiveRefinement* function, once it does not need to deal with problems related to mesh cracks avoidance in the refinement process, since all triangles will be refined.

#### Undo and Redo mechanism
If you enable the undo/redo mechanism by setting **true** to *undo* parameter in refinement functions, the following functions can be used to undo/redo a set of steps. Every function call undoes/redoes a single step.
```cpp
/**
 * Perform an undo step.
 */
void makeUndo( );

/**
 * Perform a redo step.
 */
void makeRedo( );
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
```
When you want to undo/redo all refinement steps you can do

```cpp
//Perform all undo steps.
while (!refinement->undoIsEmpty( ))
{
    s->makeUndo( );
}

...
//or
...

//Perform all redo steps.
while (!refinement->redoIsEmpty( ))
{
    s->makeRedo( );
}
```

In a real application, each undo/redo step is composed of a set of refinement operations. In this case, the application must manage an external application operations stack.
