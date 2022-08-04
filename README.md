# PyV - A toolbox for performing Particle Image Velocimetry and associated analysis

This package is geared towards collective cell migration but can also be used in other circumstances. The underlying approach is based on MATPIV[1]. In its current form the toolbox has methods for calculating the root-mean-square velocity of the field as well as the linear-order parameter (how much the vectors are aligned). There is also outputs for an alignment map (-1 negatively aligned with mean vector, 1 aligned with mean vector) and orientation. It's possible to export the fields as csv files export the rendered images for publication.

Interface:
![image](https://user-images.githubusercontent.com/45679976/182884849-68af874c-8348-48c9-a64d-fee2029fba42.png)

[1] https://www.mn.uio.no/math/english/people/aca/jks/matpiv/
