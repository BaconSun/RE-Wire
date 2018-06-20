import bpy
import bpy_extras
import numpy
from mathutils import Matrix
from mathutils import Vector

ws = "/home/bacon/Desktop/WireArt/Models/test/"

#---------------------------------------------------------------
# 3x4 P matrix from Blender camera
#---------------------------------------------------------------

# Build intrinsic camera parameters from Blender camera data
def get_calibration_matrix_K_from_blender(camd):
    f_in_mm = camd.lens
    scene = bpy.context.scene
    resolution_x_in_px = scene.render.resolution_x
    resolution_y_in_px = scene.render.resolution_y
    scale = scene.render.resolution_percentage / 100
    sensor_width_in_mm = camd.sensor_width
    sensor_height_in_mm = camd.sensor_height
    # only consider the AUTO mode. Since the pixel is square and no skew, the sensor siz should have the same
    # size of the resolution. Thus, we set the camera to have same pixles in unit length in both directions
    if (camd.sensor_fit == 'AUTO'):
        # the sensor height is fixed (sensor fit is horizontal), 
        # the sensor width is effectively changed with the pixel aspect ratio
        max_size = max(resolution_x_in_px, resolution_y_in_px)
        if max_size == resolution_x_in_px:
            s_u = resolution_x_in_px * scale / sensor_width_in_mm
            s_v = s_u
        else: 
            s_v = resolution_y_in_px * scale / sensor_width_in_mm
            s_u = s_v
    else: 
        raise Exception('Should not use mode other than Auto')


    # Parameters of intrinsic calibration matrix K
    alpha_u = f_in_mm * s_u
    alpha_v = f_in_mm * s_v
    u_0 = resolution_x_in_px * scale / 2
    v_0 = resolution_y_in_px * scale / 2
    skew = 0 # only use rectangular pixels

    K = Matrix(
        ((alpha_u, skew,    u_0),
        (    0  , alpha_v, v_0),
        (    0  , 0,        1 )))
    return K


# Returns camera rotation and translation matrices from Blender.
# 
# There are 3 coordinate systems involved:
#    1. The World coordinates: "world"
#       - right-handed
#    2. The Blender camera coordinates: "bcam"
#       - x is horizontal
#       - y is up
#       - right-handed: negative z look-at direction
#    3. The desired computer vision camera coordinates: "cv"
#       - x is horizontal
#       - y is down (to align to the actual pixel coordinates 
#         used in digital images)
#       - right-handed: positive z look-at direction
def get_3x4_RT_matrix_from_blender(cam):
    # bcam stands for blender camera
    R_bcam2cv = Matrix(
        ((1, 0,  0),
         (0, -1, 0),
         (0, 0, -1)))

    # Transpose since the rotation is object rotation, 
    # and we want coordinate rotation
    # R_world2bcam = cam.rotation_euler.to_matrix().transposed()
    # T_world2bcam = -1*R_world2bcam * location
    #
    # Use matrix_world instead to account for all constraints
    location, rotation = cam.matrix_world.decompose()[0:2]
    R_world2bcam = rotation.to_matrix().transposed()

    # Convert camera location to translation vector used in coordinate changes
    # T_world2bcam = -1*R_world2bcam*cam.location
    # Use location from matrix_world to account for constraints:     
    T_world2bcam = -1*R_world2bcam * location

    # Build the coordinate transform matrix from world to computer vision camera
    R_world2cv = R_bcam2cv*R_world2bcam
    T_world2cv = R_bcam2cv*T_world2bcam

    # put into 3x4 matrix
    RT = Matrix((
        R_world2cv[0][:] + (T_world2cv[0],),
        R_world2cv[1][:] + (T_world2cv[1],),
        R_world2cv[2][:] + (T_world2cv[2],)
         ))
    return numpy.matrix(R_world2cv), numpy.matrix(T_world2cv), RT


def get_3x4_P_matrix_from_blender(cam):
    K = get_calibration_matrix_K_from_blender(cam.data)
    R, T, RT = get_3x4_RT_matrix_from_blender(cam)
    return numpy.matrix(K*RT), numpy.matrix(K), R, T, numpy.matrix(RT)


if __name__ == "__main__":
    camera_objs = []

    for i in bpy.data.objects:
        if i.type == 'CAMERA':
            camera_objs.append(i)
            
            P, K, R, T, RT = get_3x4_P_matrix_from_blender(i)
            
            f=open(ws + i.name + ".txt",'w')
            f.close()
            f=open(ws + i.name + ".txt",'ab')
            numpy.savetxt(f, K)  # to select precision, use e.g. fmt='%.2f'
            numpy.savetxt(f, P)
            numpy.savetxt(f, R)
            numpy.savetxt(f, T)
            f.close()
            
    bpy.data.scenes['Scene'].frame_end = 2

    bpy.ops.screen.frame_jump(end=False)
    scene = bpy.data.scenes['Scene']
    render = bpy.ops.render
    scene.render.resolution_x = 640
    scene.render.resolution_y = 480
    scene.render.resolution_percentage = 100
    MAX_LEN = len(str(scene.frame_end))

    for frame in range(scene.frame_start, scene.frame_end):
        for i in camera_objs:
            scene.camera = i
            scene.render.filepath = ws + i.name + "_" + '0'*(MAX_LEN - len(str(frame))) + str(frame)
            render.render(write_still=True)
    print ("Done!")