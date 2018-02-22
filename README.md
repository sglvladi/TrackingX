# TrackingX

TrackingX is an Object Oriented MATLAB toolkit for Multi-Target Tracking, aimed at providing a common framework for swift prototyping and evaluation of multi-target tracking algorithms.  

The class architecture and interfaces of TrackingX are being designed such that algorithms of the same class (e.g. Dynamic Models) provide a common interface to all classes that utilise them (e.g. Filters), thus allowing for a virtual plug&play relation between various classes of algorithms.

Here is a temporary and partial class diagram of the intended TrackingX architecture. The full diagram should be ready soon...

Note that the letter "X" has been appended to the name of each TrackingX class in order to avoid conflicts with other matlab toolboxes. For example, classes like "KalmanFilterX" and "ParticleFilterX" are in fact implementations of the standard Kalman Filter and Particle Filter algorithms.

![Class Diagram so far...](https://github.com/sglvladi/TrackingX/blob/dev/_docs/TrackingX%20Class%20diagram%20(draft).jpg)


NOTE: Development is still in progress.. 
