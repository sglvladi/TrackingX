# TrackingX

TrackingX is an Object Oriented MATLAB toolkit for Multi-Target Tracking, aimed at providing a common framework for swift prototyping and evaluation of multi-target tracking algorithms.  

The class architecture and interfaces of TrackingX are being designed such that algorithms of the same class (e.g. Dynamic Models) provide a common interface to all classes that utilise them (e.g. Filters), thus allowing for a virtual plug&play relation between various classes of algorithms.

Below is an up-to-date class diagram of the intended TrackingX architecture.

![Class Diagram so far...](./_docs/TrackingX%20Class%20diagram%20(draft).jpg)

Note that the letter "X" has been appended to the name of each TrackingX class in order to avoid conflicts with other matlab toolboxes. For example, classes like "KalmanFilterX" and "ParticleFilterX" are in fact implementations of the standard Kalman Filter and Particle Filter algorithms.

Prerequisites
-------------
TrackingX depends on the following toolboxes:
* [Robotics System Toolbox](https://uk.mathworks.com/products/robotics.html?s_tid=AO_PR_info):
    The SystematicResamplerX and MultinomialResamplerX classes extend the matlabshared.tracking.internal.SystematicResampler and matlabshared.tracking.internal.MultinomialResampler classes, which are part of this toolbox.

NOTE: Development is still in progress.. 