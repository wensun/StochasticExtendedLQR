#pragma once

// include the types
#include "callisto41Types.h"

// the resulting manual of this file can be found at http://www.nieuwenhuisen.nl/callisto/docs/html/index.html

/*! \mainpage Callisto 4.1b
 *
 * \section intro Callisto function reference
 *
 * More information about Callisto can be found at http://www.nieuwenhuisen.nl/callisto/callisto.php
 */

/*!
  \defgroup Global Global functions
  \defgroup ObjGroup Object/group functions
    \defgroup Primitives Functions to create primitives 
    \ingroup ObjGroup
    \defgroup Animation Animation functions 
    \ingroup ObjGroup
    \defgroup Retrieval Retrieve objects/groups
    \ingroup ObjGroup
  \defgroup CollisionDetection Collision detection
  \defgroup Visualisation Visualisation
  \defgroup Statistics Statistical information
  \defgroup Errors Error return codes in Callisto
*/

/*! This function initializes Callisto, starts the output window and GUI.
	\param visualisation Set to false if you dont want a visualisation window, default is TRUE.
  \param logFile Provide a filename to create a Callisto logging file (WARNING: this may become large!).
	\returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup Global
*/
int CAL_Initialisation (bool visualisation=true, char *logFile=CAL_NULL);

/*! This function ends Callisto and cleans up memory.
	\returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup Global
*/
int CAL_End ();

/*! This function stops the visualisation until CAL_ResumeVisualisation is called.
  This function can be used to temporarily stop visualisation, for example when you need all processor power to do some calculations.
	\returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup Visualisation
*/
int CAL_SuspendVisualisation ();

/*! This function resumes the visualisation after CAL_SuspendVisualisation.
	\returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup Visualisation
*/
int CAL_ResumeVisualisation ();

/*! Shows a view. Every view has its own unique ID.
    There are at most 4 views (0..3). 0 being the main view. This view cannot be switched on/off.
  \param viewID The ID of the view.
  \param caption The caption of the window of the view, default is no caption.
  \returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup Visualisation
*/
int CAL_ShowView (int viewID, char* caption="");

/*! Hides a view. Every view has its own unique ID.
    There are at most 4 views (0..3). 0 being the main view. This view cannot be switched on/off.
  \param viewID The ID of the view.
  \returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup Visualisation
*/
int CAL_HideView (int viewID);

/*! Change the current view parameters.
  \param viewID The ID of the view.
  \param options The option you want to change. Multiple parameters can be changed at once by using the | operator. Legal values can be found \ref Global "here".
  \returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup Visualisation
*/
int CAL_SetViewOptions (int viewID, long options);

/*! Set the fog distance. This value is inpendent of the view.
  \param dist The new fog distance.
  \returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup Visualisation
*/
int CAL_SetViewFogDistance (CAL_scalar dist);

/*! Set the navigation speed. This is the speed with the keyboard navigation and middle mouse button work.
  \param viewID The ID of the view.
  \param dist The new navigation speed.
  \returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup Visualisation
*/
int CAL_SetViewNavigationSpeed (int viewID, CAL_scalar dist);

/*! Set the near clipping distance.
  \param viewID The ID of the view.
  \param dist The new near clipping distance, legal values are between 0 and 100.
  \returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup Visualisation
*/
int CAL_SetViewNearClippingDistance (int viewID, CAL_scalar dist);

/*! This function sets the position and orientation of the camera in the visualisation window.
    It is positioned at eye and looks to the point defined with view. By default the up vector is (0,1,0)
    You can change this if necessary.
  \param viewID The ID of the view to set the view parameters for, set to 0 if you don't know what this is.
  \param eyeX The x position of the camera.
  \param eyeY The y position of the camera.
  \param eyeZ The z position of the camera.
  \param lookAtX The x position of the look at point.
  \param lookAtY The y position of the look at point.
  \param lookAtZ The z position of the look at point.
  \param roll The roll about the eye-view vector.
	\returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup Visualisation
*/
int CAL_SetViewParams (int viewID, CAL_scalar eyeX, CAL_scalar eyeY, CAL_scalar eyeZ, CAL_scalar lookAtX, CAL_scalar lookAtY, CAL_scalar lookAtZ, CAL_scalar roll=0);

/*! This function gets the current position and orientation of the camera in a visualisation window.
    It is positioned at eye and looks to the point defined with view. 
  \param viewID The ID of the view get the view parameters from.
  \param *eyeX The x position of the camera.
  \param *eyeY The y position of the camera.
  \param *eyeZ The z position of the camera.
  \param *lookAtX The x position of the look at point.
  \param *lookAtY The y position of the look at point.
  \param *lookAtZ The z position of the look at point.
  \param *roll The roll about the eye-view vector.
	\returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup Visualisation
*/
int CAL_GetViewParams (int viewID, CAL_scalar *eyeX, CAL_scalar *eyeY, CAL_scalar *eyeZ, CAL_scalar *lookAtX, CAL_scalar *lookAtY, CAL_scalar *lookAtZ, CAL_scalar *roll);

/*! Adds a keystate for the view animation.
  \param viewID The ID of the view.
  \param time The time stamp for which the configuration accounts.
  \param eye Three scalar values that determine the position of the eye point.
  \param view Three scalar values that determine the position of the look-at point.
  \param roll The roll about the eye-view vector.
	\returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup Animation
*/
int CAL_AddViewKeyState (int viewID, CAL_scalar time, CAL_scalar *eye=CAL_NULL, CAL_scalar *view=CAL_NULL, CAL_scalar roll=-1);

/*! Determines whether a view animation is cyclic or not.
  \param viewID The ID of the view.
  \param cyclic Set to true for a cyclic motion, or false for a non-cyclic motion.
	\returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup Animation
*/
int CAL_SetViewCyclic (int viewID, bool cyclic);

/*! Deletes all existing view key states.
  \param viewID The ID of the view.
	\returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup Animation
*/
int CAL_ClearViewKeyStates (int viewID);

/*! Sets the background color of the visualisation window.
  \param viewID The ID of the view get the view parameters from, set to 0 if you don't know what this is.
  \param red The red component of the color (0...1).
  \param green The green component of the color (0...1).
  \param blue The blue component of the color (0...1).
	\returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup Visualisation
*/
int CAL_SetBackgroundColor (int viewID, CAL_scalar red, CAL_scalar green, CAL_scalar blue);

/*! Adds a resource (i.e. a directory) to the list of available resources.
  \param pathName The name of the path to add to the list of resources.
  \ingroup Global
*/
int CAL_AddTextureResource (char* pathName);

/*! Loads a texture to memory. The texture must be in the .ppm or the .png format. Do not forget to add a texture resource path first by using CAL_AddTextureResource.
  \param textureID The ID you want to give the texture, there is room for 500 textures numbered (0..499).
  \param fileName The filename of the texture. Do not add the filepath (this should be added through CAL_AddTextureResource). The texture has to be in .ppm, .png or .jpg format.
  \returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup Global
*/
int CAL_LoadTexture (int textureID, char* fileName);

/*! Add a texture from memory.
  \param textureID The ID you want to give the texture.
  \param width The width in pixels of the texture. Must be a power of 2.
  \param height The height in pixels of the texture. Must be a power of 2.
  \param tex The textures in RGB format. The format is 8 bits red, 8 bits green and 8 bits blue (3 bytes per pixel). Length is width*height*3.
  \returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup Global
*/
int CAL_SetTextureFromMem (int textureID, int width, int height, unsigned char *tex);

/*! Saves the content of a view in .bmp format to disk.
  \param viewID The ID of the view to save, set to 0 if you don't know what this is.
  \param fileName The file name of the .bmp.
  \returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup Visualisation
*/
int CAL_ScreenCapture (int viewID, char* fileName);

/*! Load a scene in XML or VRML format from disk.
  \param *fileName The name of the file to load.
  \param parentID The parent group to put the loaded scene in, use 0 for no parent.
  \param **error String with possible error string. Can be omitted.
	\returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup Global
*/
int CAL_LoadScene (char* fileName, int parentID, char **error=CAL_NULL);

/*! Save a (part of a) scene in XML or VRML format from disk.
  \param *fileName The name of the file to load.
  \param viewID The view id which needs to be saved (only visible objects are saved)
  \param groupID The group id of the group that needs to be saved
	\returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup Global
*/
int CAL_SaveScene (char* fileName, int viewID=0, int groupID=0);

/*! Change the time, all dynamic groups and objects will adapt to this time.
  \param time The time;
	\returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup Animation
*/
int CAL_SetTime (CAL_scalar time);

/*! Creates an object group.
  \param groupID This is set to the group ID the group gets.
  \param parentID The parent group of this group, use 0 if this group does not need a parent.
  \param colCheck True if this group is used for collision checks, false if not.
  \param name The name of the group as shown in the GUI, default is no name.
  \param collapsed The group will appear collapsed in the interface.
  \param visible Set to either CAL_TRUE, CAL_FALSE or CAL_USEPARENT.
	\returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup ObjGroup
*/
int CAL_CreateGroup (int* groupID, int parentID, bool colCheck, char* name="", bool collapsed=false, int visible=CAL_USEPARENT);

/*! Delets a group, its child groups and all its objects.
  \param groupID The ID of the group to delete.
	\returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup ObjGroup
*/
int CAL_DestroyGroup (int groupID);

/*! Deletes a groups child groups and all its objects.
  \param groupID The ID of the group to empty.
  \param subGroups Set to TRUE if subgroups of group should also be removed, default is FALSE.
	\returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup ObjGroup
*/
int CAL_EmptyGroup (int groupID, bool subGroups=false);

/*! Move a group to a new parent.
  \param groupID The ID of the group to move.
  \param parentID The ID of the new parent.
	\returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup ObjGroup
*/
int CAL_MoveGroup (int groupID, int parentID);

/*! Place a group and all its objects at a new position.
  \param groupID The ID of the group to translate.
  \param x The x component of the new position.
  \param y The y component of the new position.
  \param z The z component of the new position.
	\returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup ObjGroup
*/
int CAL_SetGroupPosition (int groupID, CAL_scalar x, CAL_scalar y, CAL_scalar z);

/*! Place a group and all its objects in a new orientation using Euler angles (in radians).
  \param groupID The ID of the group to translate.
  \param xRot The rotation orientation with respect to the x-axis;
  \param yRot The rotation orientation with respect to the y-axis;
  \param zRot The rotation orientation with respect to the z-axis;
	\returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup ObjGroup
*/
int CAL_SetGroupOrientation (int groupID, CAL_scalar xRot, CAL_scalar yRot, CAL_scalar zRot);

/*! Scale a group and all its objects.
  \param groupID The ID of the group to scale.
  \param xScale The scaling factor in the x direction.
  \param yScale The scaling factor in the y direction.
  \param zScale The scaling in the z direction.
	\returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup ObjGroup
*/
int CAL_SetGroupScaling (int groupID, CAL_scalar xScale, CAL_scalar yScale, CAL_scalar zScale);

/*! Rotate a group and all its objects using a quaternion.
  \param groupID The ID of the group to rotate.
  \param x The x component of the quaternion.
  \param y The y component of the quaternion.
  \param z The z component of the quaternion.
  \param w The w component of the quaternion.
	\returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup ObjGroup
*/
int CAL_SetGroupQuaternion (int groupID, CAL_scalar x, CAL_scalar y, CAL_scalar z, CAL_scalar w);

/*! Set the orientation matrix of a group.
  \param groupID The ID of the group to rotate.
  \param mat The 3x3 matrix to set the orientation.
	\returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup ObjGroup
*/
int CAL_SetGroupMatrix (int groupID, CAL_matrix3 mat);

/*! Spherically expand a group. Collision will occur with the set of points whose distance to the objects is at most the clearance.
    This function overrules individual group clearance.
  \param groupID The ID of the group to translate.
  \param c The clearance.
	\returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup ObjGroup
*/
int CAL_SetGroupClearance (int groupID, CAL_scalar c);

/*! Sets the color of a group.
  \param groupID The ID of the group.
  \param red The red component of the color (0...1).
  \param green The green component of the color (0...1).
  \param blue The blue component of the color (0...1).
  \param alpha The alpha value of the color (0...1), 0 is fully transparant, 1 is fully opaque. Default is 1.
  \param receiveShadows Set when the objects need to receive shadows from other objects.
  \param subGroups Set to TRUE if subgroups of group should also get new color, default is FALSE.
  \param sID The ID of the material which gets this color, default is material 0.
	\returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup ObjGroup
*/
int CAL_SetGroupColor (int groupID, CAL_scalar red, CAL_scalar green, CAL_scalar blue, CAL_scalar alpha=1, bool receiveShadows=true, bool subGroups=false, int sID=0);

/*! Set the texture for a group.
  \param groupID The ID of the group.
  \param textureID The ID of the texture.
  \param xtile The number of times to repeat the texture in the x-direction.
  \param ytile The number of times to repeat the texture in the y-direction.
  \param alpha The alpha value of the color (0...1), 0 is fully transparant, 1 is fully opaque. Default is 1.
  \param receiveShadows Set when the objects need to receive shadows from other objects.
  \param subGroups Set to TRUE if subgroups of group should also get new color. Default is FALSE.
  \param sID The ID of the material which gets this texture, default is material 0.
  \returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup ObjGroup
*/
int CAL_SetGroupTexture (int groupID, int textureID, CAL_scalar xtile, CAL_scalar ytile, CAL_scalar alpha=1, bool receiveShadows=true, bool subGroups=false, int sID=0);

/*! Set whether the objects of the group cast shadows. The default value is true.
    Casting is shadow is independent of the color/texture the object has.
  \param groupID The ID of the group.
  \param castShadows Boolean the states whether the objects should cast shadows or not.
  \param subGroups Set to TRUE if casting shadows of subgroups should also be set. Default is FALSE.
  \returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup ObjGroup
*/
int CAL_SetGroupCastShadows (int groupID, bool castShadows, bool subGroups=false);

/*! Sets the active material of the group.
  \param groupID The ID of the group
  \param subGroups Set to true if subgroups also have to change the active material.
  \param sID The ID of the new active material.
  \ingroup ObjGroup
*/
int CAL_SetGroupActiveMaterial (int groupID, bool subGroups, int sID);

/*! Sets the collision check capability of a group.
  \param groupID The ID of the group.
  \param colCapable True if the group can be collision checked.
  \param subGroups Set to true if subgroups also have to change their collision check capability recursively.
	\returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup ObjGroup
*/
int CAL_SetGroupCollisionCheckCapability (int groupID, bool colCapable, bool subGroups=false);

/*! Sets the visibility of a group in a particular view. This does not effect collision checks.
  \param groupID The ID of the group.
  \param viewID The ID of the view.
  \param visible True if the group should be visible. False for invisibility.
  \param subGroups Set to true if subgroups also have to change their visibility recursively.
	\returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup ObjGroup
*/
int CAL_SetGroupVisibility (int groupID, int viewID, bool visible, bool subGroups=false);

/*! Change the name of a group.
  \param groupID The ID of the group.
  \param name The new name
  \ingroup ObjGroup
*/
int CAL_SetGroupName (int groupID, char *name);

/*! Adds a keystate for the group animation.
  \param groupID The ID of the group.
  \param time The time stamp for which the configuration accounts.
  \param position Three scalar values that determine the position of the group at this keystate.
  \param orientation The orientation of the group, this can consist of Euler angles or a quaternion, depending on the last parameter.
  \param scaling The scaling of the group, consisting of three values.
  \param isQuat Set this to false if the orientation parameter consists of Euler angles.
  \param visible Can be either CAL_TRUE, CAL_FALSE or CAL_NOTSET. If not set, the group value (by CAL_SetGroupVisibility) is used for the whole animation.
	\returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup Animation
*/
int CAL_AddGroupKeyState (int groupID, CAL_scalar time, CAL_scalar *position, CAL_scalar *orientation=CAL_NULL, CAL_scalar *scaling=CAL_NULL, bool isQuat=true, int visible=CAL_NOTSET);

/*! Change the group motion parameters.
  \param groupID The ID of the group.
  \param options Legal values are: CAL_CYCLIC, CAL_NONCYCLIC, CAL_INTERPOLATION, CAL_NOINTERPOLATION.
  \returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup Animation
*/
int CAL_SetGroupMotionOptions (int groupID, int options);

/*! Deletes all existing group key states.
  \param groupID The ID of the group.
  \param subGroups Set to TRUE if subgroups of group should also be cleared.
	\returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup Animation
*/
int CAL_ClearGroupKeyStates (int groupID, bool subGroups=false);

/*! Clones a group including all objects.
  \param groupIDNew Set to the group ID of the clone.
  \param groupID The ID of the group to clone.
  \param parentID The parent group of the clone (0 for no parent).
  \param colCheck Set to false if clone group does not need collision checks.
  \param name The name of the group, default is no name.
  \param cloneObjs Optional flag that indicates whether to clone the objects in the group as well.
  \returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup ObjGroup
*/
int CAL_CloneGroup (int* groupIDNew, int groupID, int parentID, bool colCheck, char* name="", bool cloneObjs=true);

/*! Clones a group including all objects and child groups.
  \param newgroupIDs List of new group IDs, can be NULL.
  \param groupID The ID of the group to clone.
  \param parentID The parent group of the clone (0 for no parent).
  \param nr The size of the ids and names lists, can be 0.
  \param ids List of id's corresponding to the names list, can be NULL.
  \param names Optional list of new names for the clone, defaults to NULL.
  \param cloneObjs Optional flag that indicates whether to clone the objects in the group as well, defaults to true.
  \param keepColCap Optional flag that indicates whether the collision capabilities should be preserved. Default clones do not have collision check capabilities.
	\returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup ObjGroup
*/
int CAL_CloneGroupRecursive (int* newgroupIDs, int groupID, int parentID, int nr, int* ids, char** names=CAL_NULL, bool cloneObjs=true, bool keepColCap=false);

/*! Move an object to another group.
  \param objID The ID of the object to destroy.
  \param groupID The ID of the new group.
	\returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup ObjGroup
*/
int CAL_MoveObject (int objID, int groupID);

/*! Destroys an object.
  \param objID The ID of the object to destroy.
	\returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup ObjGroup
*/
int CAL_DestroyObject (int objID);

/*! Place an object at a new position.
  \param objID The ID of the object to translate.
  \param x The x component of the new position.
  \param y The y component of the new position.
  \param z The z component of the new position.
	\returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup ObjGroup
*/
int CAL_SetObjectPosition (int objID, CAL_scalar x, CAL_scalar y, CAL_scalar z);

/*! Place an object in a new orientation using Euler angles (in radians).
  \param objID The ID of the object to translate.
  \param xRot The orientation with respect to the x-axis;
  \param yRot The orientation with respect to the y-axis;
  \param zRot The orientation with respect to the z-axis;
	\returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup ObjGroup
*/
int CAL_SetObjectOrientation (int objID, CAL_scalar xRot, CAL_scalar yRot, CAL_scalar zRot);

/*! Scale an object.
  \param objID The ID of the object to scale.
  \param xScale The scaling factor in the x direction.
  \param yScale The scaling factor in the y direction.
  \param zScale The scaling factor in the z direction.
	\returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup ObjGroup
*/
int CAL_SetObjectScaling (int objID, CAL_scalar xScale, CAL_scalar yScale, CAL_scalar zScale);

/*! Rotate an object using a quaternion.
  \param objID The ID of the object to rotate.
  \param x The x component of the quaternion.
  \param y The y component of the quaternion.
  \param z The z component of the quaternion.
  \param w The w component of the quaternion.
	\returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup ObjGroup
*/
int CAL_SetObjectQuaternion (int objID, CAL_scalar x, CAL_scalar y, CAL_scalar z, CAL_scalar w);

/*! Set the orientation matrix of an object.
  \param objID The ID of the object to rotate.
  \param mat The 3x3 matrix to set the orientation.
	\returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup ObjGroup
*/
int CAL_SetObjectMatrix (int objID, CAL_matrix3 mat);

/*! Set the WORLD matrix of an object. This works even when group matrix is set.
  \param objID The ID of the object.
  \param matrix Should be of type CAL_matrix;
  \returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup ObjGroup
*/
int CAL_SetObjectWorldMatrix (int objID, CAL_matrix4 *matrix);

/*! Spherically expand an object. Collision will occur with the set of points whose distance to the object is at most the clearance.
    If CAL_SetGroupClearance is called later, individual object clearance will be overruled.
  \param objID The ID of the object to translate.
  \param clearance The clearance.
	\returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup ObjGroup
*/
int CAL_SetObjectClearance (int objID, CAL_scalar clearance);

/*! Sets the color of an object.
  \param objID The ID of the object.
  \param red The red component of the color (0...1).
  \param green The green component of the color (0...1).
  \param blue The blue component of the color (0...1).
  \param alpha The alpha value of the color (0...1), 0 is fully transparant, 1 is fully opaque. Default is 1.
  \param receiveShadows Set when the object needs to receive shadows from other objects.
  \param sID The ID of the material which gets this color. Default value is material 0.
	\returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup ObjGroup
*/
int CAL_SetObjectColor (int objID, CAL_scalar red, CAL_scalar green, CAL_scalar blue, CAL_scalar alpha=1, bool receiveShadows=true, int sID=0);

/*! Set the texture for an object.
  \param objID The ID of the object.
  \param textureID The ID of the texture.
  \param xtile The number of times to repeat the texture in the x-direction.
  \param ytile The number of times to repeat the texture in the y-direction.
  \param alpha The alpha value of the color (0...1), 0 is fully transparant, 1 is fully opaque. Default is 1.
  \param receiveShadows Set when the object needs to receive shadows from other objects.
  \param sID The ID of the material which gets this texture, default is material 0.
  \returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup ObjGroup
*/
int CAL_SetObjectTexture (int objID, int textureID, CAL_scalar xtile, CAL_scalar ytile, CAL_scalar alpha=1, bool receiveShadows=true, int sID=0);

/*! Set whether the object casts shadow. The default value is true.
    Casting is shadow is independent of the color/texture the object has.
  \param objID The ID of the object.
  \param castShadows Boolean the states whether the object should cast shadows or not.
  \returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup ObjGroup
*/
int CAL_SetObjectCastShadows (int objID, bool castShadows);

/*! Sets the active material of the group.
  \param objID The ID of the object.
  \param sID The ID of the new active material.
  \ingroup ObjGroup
*/
int CAL_SetObjectActiveMaterial (int objID, int sID);

/*! Adds a keystate for the object animation.
  \param objID The ID of the object.
  \param time The time stamp for which the configuration accounts.
  \param position Three scalar values that determine the position of the object at this keystate.
  \param orientation The orientation of the object, this can consist of Euler angles or a quaternion, depending on the last parameter.
  \param scaling The scaling of the group, consisting of three values.
  \param isQuat Set this to false if the orientation parameter consists of Euler angles.
	\returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup Animation
*/
int CAL_AddObjectKeyState (int objID, CAL_scalar time, CAL_scalar *position, CAL_scalar *orientation=CAL_NULL, CAL_scalar *scaling=CAL_NULL, bool isQuat=true);

/*! Change the objectmotion parameters.
  \param objID The ID of the object.
  \param options Legal values are: CAL_CYCLIC, CAL_NONCYCLIC, CAL_INTERPOLATION, CAL_NOINTERPOLATION.
  \returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup Animation
*/
int CAL_SetObjectMotionOptions (int objID, int options);

/*! Deletes all existing object key states.
  \param objID The ID of the object.
	\returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup Animation
*/
int CAL_ClearObjectKeyStates (int objID);

/*! Returnes the ID of a group or object with a certain name.
  \param ID This is set to the ID.
  \param name The name of the group/object.
	\returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup Retrieval
*/
int CAL_GetID (int* ID, char* name);

/*! Get the properties of a group.
  \param groupID The ID of the group.
  \param CALGroup Pointer to an SCALGroup-structure.
  \returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup Retrieval
*/
int CAL_GetGroup (int groupID, void *CALGroup);

/*! Returnes the ID of the nr'th childgroup.
  \param groupID The ID of the group.
  \param nr The nr of the childgroup.
  \param childGroupID This value is set to the nr'th childgroup.
	\returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup Retrieval
*/
int CAL_GetGroupChildID (int groupID, int nr, int *childGroupID);

/*! Returnes the ID of the nr'th object.
  \param groupID The ID of the group.
  \param nr The nr of the object.
  \param objectID This value is set to the nr'th object.
	\returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup Retrieval
*/
int CAL_GetGroupObjectID (int groupID, int nr, int *objectID);


/*! Get the type of an object (#CAL_BOX, #CAL_CYLINDER etc.).
  \param objID The ID of the object.
  \param objType This will be set to the object type according to the values in callistoTypes.h.
  \returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup Retrieval
*/
int CAL_GetObjectType (int objID, int* objType);

/*! Get the WORLD matrix of an object.
  \param objID The ID of the object.
  \param matrix Should be of type CAL_matrix;
  \returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup Retrieval
*/
int CAL_GetObjectWorldMatrix (int objID, CAL_matrix4 *matrix);

/*! Get the properties of an object.
  \param objID The ID of the object.
  \param SCALObj Pointer to an SCAL-object (SCALBox, SCALSphere etc.). This has to be of the right type (the type can be retrieved by using CAL_GetObjectType).
         Note that to retrieve a CAL_ELEVATIONGRID object, SCALPolygonGroup needs to be used.
  \returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup Retrieval
*/
int CAL_GetObject (int objID, void* SCALObj);

/*! Set a callback function, this is called when the user pressed a key in a view.
  \param cb The adress of the callback function.
	\returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup Global
*/
int CAL_SetKeypressCallback (CAL_KeypressCallback cb);

/*! Set a callback function, this is called when the user selects an object by clicking the left mouse button while pressing SHIFT.
  \param cb The adress of the callback function.
	\returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup Global
*/
int CAL_SetObjectSelectCallback (CAL_ObjectSelectCallback cb);

/*! Create a box.
  \param groupID The group ID to put the object in.
  \param xw The width of the object.
  \param yw The height of the object.
  \param zw The depth of the object.
  \param x The x position of the object.
  \param y The y position of the object.
  \param z The z position of the object.
  \param objID Set to the object ID.
  \param name Name of the object. Default is no name.
	\returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup Primitives
*/
int CAL_CreateBox (int groupID, CAL_scalar xw, CAL_scalar yw, CAL_scalar zw, CAL_scalar x, CAL_scalar y, CAL_scalar z, int* objID=CAL_NULL, char* name="");

/*! Create a sphere.
  \param groupID The group ID to put the object in.
  \param radius The radius of the object.
  \param x The x position of the object.
  \param y The y position of the object.
  \param z The z position of the object.
  \param objID Set to the object ID.
  \param name Name of the object. Default is no name.
	\returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup Primitives
*/
int CAL_CreateSphere (int groupID, CAL_scalar radius, CAL_scalar x, CAL_scalar y, CAL_scalar z, int* objID=CAL_NULL, char* name="");

/*! Create a cylinder.
  \param groupID The group ID to put the object in.
  \param radius The radius of the object.
  \param height The height of the object.
  \param x The x position of the object.
  \param y The y position of the object.
  \param z The z position of the object.
  \param objID Set to the object ID.
  \param name Name of the object. Default is no name.
	\returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup Primitives
*/
int CAL_CreateCylinder (int groupID, CAL_scalar radius, CAL_scalar height, CAL_scalar x, CAL_scalar y, CAL_scalar z, int* objID=CAL_NULL, char* name="");

/*! Create a sphere.
  \param groupID The group ID to put the object in.
  \param radius The radius of the bottom of the object.
  \param height The height of the object.
  \param x The x position of the object.
  \param y The y position of the object.
  \param z The z position of the object.
  \param objID Set to the object ID.
  \param name Name of the object. Default is no name.
	\returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup Primitives
*/
int CAL_CreateCone (int groupID, CAL_scalar radius, CAL_scalar height, CAL_scalar x, CAL_scalar y, CAL_scalar z, int* objID=CAL_NULL, char* name="");

/*! Create a triangle mesh.
  \param groupID The group ID to put the object in.
  \param nrTriangles The number of triangles the object consists of.
  \param *p List of coordinates of the object. The size of the list must be 9*nrTriangles.
  \param *texCoords List of texture coordinates. Each texture coordinate consists of two values betwee 0 and1. 
          There should be 2 coordinates for each point in the list. Thus its size should be 6*nrTriangles.
  \param objID Set to the object ID.
  \param name Name of the object. Default is no name.
	\returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup Primitives
*/
int CAL_CreateTriangles (int groupID, int nrTriangles, CAL_scalar* p, CAL_scalar* texCoords=CAL_NULL, int* objID=CAL_NULL, char* name="");

/*! Create a polyline. This is for the visualisation only.
  \param groupID The group ID to put the object in.
  \param nl The number of lines.
  \param *np The number of points each line consists of.
  \param *p List of coordinates of the lines. Its size should be is 3 * np * nl.
  \param objID Set to the object ID.
  \param name Name of the object. Default is no name.
	\returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup Primitives
*/
int CAL_CreatePolyline (int groupID, int nl, int *np, CAL_scalar *p, int* objID=CAL_NULL, char* name="");

/*! Create a tetrahedron.
  \param groupID The group ID to put the object in.
  \param *p List of coordinates of the object. The size of the list must be 3*4 since a tetrahydron consists of 4 points. 
    The first 3 points define the ground plane, the fourth is the top. The ground plane should be defined in a counter-clockwise order.
  \param objID Set to the object ID.
  \param name Name of the object. Default is no name.
	\returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup Primitives
*/
int CAL_CreateTetrahedron (int groupID, CAL_scalar* p, int* objID=CAL_NULL, char* name="");

/*! Create an elevation grid on the XZ plane.
  \param groupID The group ID to put the object in.
  \param xDim The number of x coordinates.
  \param zDim The number of z coordinates.
  \param xStep The stepsize along the x-axis.
  \param zStep The stepsize along the z-axis.
  \param heights The height parameters. The number of elements should be (xDim+1)*(zDim+1).
  \param objID Set to the object ID.
  \param name Name of the object. Default is no name.
	\returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup Primitives
*/
int CAL_CreateElevationGrid (int groupID, int xDim, int zDim, CAL_scalar xStep, CAL_scalar zStep, CAL_scalar* heights, int *objID=CAL_NULL, char* name="");

/*! Create a label attached to an object. The label consists of a line and a text. The text is always faced towards the viewer.
  \param objID The ID of the object.
  \param caption The caption of the label (use \n to start a new line).
  \param charHeight Size of the label text (between 0 and 100).
  \param x The x position of the label w.r.t. the object.
  \param y The y position of the label w.r.t. the object.
  \param z The z position of the label w.r.t. the object.
  \returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup ObjectDef
*/
int CAL_CreateLabel (int objID, char *caption, CAL_scalar charHeight, CAL_scalar x, CAL_scalar y, CAL_scalar z);

/*! Change the caption of the label.
  \param objID The ID of the object.
  \param caption The caption of the label.
  \returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup ObjectDef
*/
int CAL_SetLabelCaption (int objID, char *caption);

/*! Change the color of the label.
  \param objID The ID of the object.
  \param red The red component of the color.
  \param green The green component of the color.
  \param blue The blue component of the color.
  \param alpha The alpha transparency component of the color.
  \returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup ObjectDef
*/
int CAL_SetLabelColor (int objID, CAL_scalar red, CAL_scalar green, CAL_scalar blue, CAL_scalar alpha);

/*! Destroy a label.
  \param objID The ID of the object which label should be desroyed.
  \returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup ObjectDef
*/
int CAL_DestroyLabel (int objID);

/*! Create an overlay. An overlay is a box that appears at at a static position in the visualisation window. It can be used to display user defined information.
    An example of an overlay is the status overlay.
    \param overlayID The ID of the overlay.
    \param viewID The ID of the view the overlay should appear in.
    \param xPos The x coordinate of the screen position the overlay should appear.
    \param yPos The y coordinate of the screen position the overlay should appear.
    \param xDim The width of the overlay.
    \param yDim The height of the overlay.
    \param referencePos The reference from which xPos and yPos are calculated. The value can be CAL_TOPLEFT, CAL_TOPRIGHT, CAL_BOTTOMRIGHT or CAL_BOTTOMLEFT.
    \param nrLines The number of lines the overlay consists of.
    \ingroup Visualisation
*/
int CAL_CreateOverlay (int *overlayID, int viewID, int xPos, int yPos, int xDim, int yDim, int referencePos=CAL_TOPLEFT, int nrLines=1);

/*! Destroy an overlay.
  \param overlayID The ID of the overlay that needs to be destroyed.
  \ingroup Visualisation
*/
int CAL_DestroyOverlay (int overlayID);

/*! Change the overlay text.
  \param overlayID The ID of the overlay.
  \param text An array of text lines. The length of the array should be equal as the number of lines used in CAL_CreateOverlay.
  \ingroup Visualisation
*/
int CAL_SetOverlayText (int overlayID, char **text);

/*! Change the visibility of an overlay.
  \param overlayID The ID of the overlay.
  \param visibility Boolean indicating visibility.
  \ingroup Visualisation
*/
int CAL_SetOverlayVisibility (int overlayID, bool visibility);

/*! Check whether a point collides with a group.
  \param groupID The group ID to check with.
  \param x The x coordinate of the point to check.
  \param y The y coordinate of the point to check.
  \param z The z coordinate of the point to check.
  \param multiple Flag whether to find all results, or just the first (faster).
  \param *nrCols The number of collisions.
	\returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup CollisionDetection
*/
int CAL_CheckPointCollision (int groupID, CAL_scalar x, CAL_scalar y, CAL_scalar z, bool multiple, int *nrCols);

/*! Check whether a line collides with a group.
  \param groupID The group ID to check with.
  \param x0 The x coordinate of the first point to check.
  \param y0 The y coordinate of the first point to check.
  \param z0 The z coordinate of the first point to check.
  \param x1 The x coordinate of the second point to check.
  \param y1 The y coordinate of the second point to check.
  \param z1 The z coordinate of the second point to check.
  \param multiple Flag whether to find all results, or just the first (faster).
  \param *nrCols The number of collisions.
	\returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup CollisionDetection
*/
int CAL_CheckLineCollision (int groupID, CAL_scalar x0, CAL_scalar y0, CAL_scalar z0, CAL_scalar x1, CAL_scalar y1, CAL_scalar z1, bool multiple, int *nrCols);

/*! Check whether two groups collide. Groups cannot be each others subgroups.
  \param group0 The ID of the first group.
  \param group1 The ID of the second group.
  \param multiple Flag whether to find all results, or just the first (faster).
  \param *nrCols The number of collisions.
	\returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup CollisionDetection
*/
int CAL_CheckGroupCollision (int group0, int group1, bool multiple, int *nrCols);

/*! Get the two positions where the distance between two groups is smallest.
  \param groupID0 The ID of the first group.
  \param groupID1 The ID of the second group.
  \param nrPairs The number of results.
  \returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup CollisionDetection
*/
int CAL_GetClosestPairs (int groupID0, int groupID1, int *nrPairs);

/*! Get the two positions where the penetration of two groups is largest.
  \param groupID0 The ID of the first group.
  \param groupID1 The ID of the second group.
  \param nrPairs The number of results.
  \returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup CollisionDetection
*/
int CAL_GetPenetrationDepths (int groupID0, int groupID1, int *nrPairs);

/*! Get the results of the objects involved in the last collision check/penetration depth/closest pair.
  \param userResults A list of type SCALResult which contains the results. The client is reponsible of creating the list with the right size (size == count). The results for closest pairs and penetration depths are sorted on distance.
  \returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup CollisionDetection
*/
int CAL_GetResults (void* userResults);

/*! Enable statistics. Statistical information will be gathered for collision checks, closest pairs and penetration depth. Both the number of calls and the total time spent in these functions will be administered. Beware that to gather the statistics themselved a small amount of time is spent.
  \param enable This value has to be either equal to either CAL_ENABLESTATISTICS or CAL_DISABLESTATISTICS.
  \ingroup Statistics
*/
int CAL_GatherStatistics (int enable);

/*! Reset all statistics to 0.
  \ingroup Statistics
*/
int CAL_ResetStatistics ();

/*! Save the statistical information to a comma delimited file.
  \param groupID The ID of the group.
  \param fileName The path and name of the file to write the information to.
  \ingroup Statistics
*/
int CAL_SaveGroupStatistics (int groupID, char *fileName);

/*! \file callisto41.h
\brief All functions and constants of Callisto.
*/