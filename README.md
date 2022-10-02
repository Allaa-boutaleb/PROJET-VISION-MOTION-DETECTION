# PROJET-VISION-MOTION-DETECTION
A motion detection project made in C with an emphasis on low-level optimization. This project was made during my 3rd year of Computer Science at Sorbonne University, in collaboration with NASSIM AHMED ALI and under the supervision of QUENTIN MEUNIER.

Detecting moving objects in a sequence of images is a important low-level task for many applications ofcomputer vision such as video surveillance, traffic monitoring, traffic recognition, signaling language…etc. When the camera is stationary, a class of methods that is usually employed is background noise subtraction. the principle of these methods is to build a model of the scene static (i.e. without moving objects) called the background, then compare each image of the sequence to this background in order to map out moving regions, this part is called the foreground.

As part of this project, our mission was to create an image processing chain to isolate moving pixels,
to optimize this chain as much as possible via various methods elaborated in the PDF report and finally to compare these different methods when it comes to speed and memory occupation.
