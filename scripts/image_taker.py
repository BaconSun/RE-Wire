#!/usr/bin/python

import rospy
from sensor_msgs.msg import Image
import cv2
import sys
from cv_bridge import CvBridge, CvBridgeError


bridge = CvBridge()


def callback(msg):
    try:
        cv_image = bridge.imgmsg_to_cv2(msg, "mono8")
        cv2.imwrite('/home/bacon/catkin_ws/src/rewire/scripts/' +
                    sys.argv[1] + '_' + sys.argv[2] + '.jpg', cv_image)
    except CvBridgeError as e:
        print(e)


def image_taker():
    rospy.init_node('test_image_pub', anonymous=True)
    rospy.Subscriber(('/camera/image/rgb_612205000197' if sys.argv[2] == '0' else '/camera/image/rgb_610201004117'),
                     Image, callback)
    rospy.spin()


if __name__ == '__main__':
    try:
        image_taker()
    except rospy.ROSInterruptException:
        pass
