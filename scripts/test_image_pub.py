#!/usr/bin/python

import rospy
from sensor_msgs.msg import Image
import cv2
from cv_bridge import CvBridge, CvBridgeError


bridge = CvBridge()
cv_image = cv2.imread('/home/bacon/catkin_ws/src/rewire/scripts/view1_bd0.jpg', 0)
msg = bridge.cv2_to_imgmsg(cv_image, "mono8")


def talker():
    pub = rospy.Publisher('/test_image', Image, queue_size=10)
    rospy.init_node('test_image_pub', anonymous=True)
    rate = rospy.Rate(1)  # 10hz
    while not rospy.is_shutdown():
        try:
            pub.publish(msg)
        except CvBridgeError as e:
            print(e)
        rate.sleep()


if __name__ == '__main__':
    try:
        talker()
    except rospy.ROSInterruptException:
        pass
