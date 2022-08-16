# PyV - A toolbox for performing Particle Image Velocimetry and associated analysis

This package is geared towards collective cell migration but can also be used in other circumstances. The underlying approach for calculating the PIV fields is based on MATPIV[1] and is relatively fast through multiprocessing. In its current form the toolbox has methods for calculating the root-mean-square velocity of the field as well as the linear-order parameter (how much the vectors are aligned). There is also outputs for an alignment map (-1 negatively aligned with mean vector, 1 aligned with mean vector) and orientation. There is also limited support at the moment for calculating the rotational components about a point in the field.

It's possible to export the fields as either csv or hdf5 files and export the rendered images for publication.

Interface:
<img width="1435" alt="image" src="https://user-images.githubusercontent.com/45679976/184831031-df05b11b-5726-4ea1-8652-6b3678245c3c.png">

[1] https://www.mn.uio.no/math/english/people/aca/jks/matpiv/
