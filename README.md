# LoGConv
Applies the Laplacian-of-Gaussian edge-detection filter to pictures in various image editors.

This implementation started in 1998 while working at my thesis; after 20 years I ported it to:
* Gimp 2.x
* paint.net

## The Laplacian-of-Gaussian filter

The logconv filter implements an algorithm for edge detection developed in 1989 by Sotak and Boyer (read the article in Computer Vision, Graphics, and Image Processing 48, pgg. 147-189).

The Laplacian-of-Gaussian (LoG) operator has been suggested by Marr and Hildreth whilst studying the human visual system.
I will not explain you all the mathematics behind the development of this interesting filter (which you can read in the article cited above).
I will instead present you a few example of the use of this new operator and introduce you to the correct use of the parameters required by the filter.

First of all, a perspective plot of the LoG operator.

![loggraph](/docs/loggraph.jpeg)

The Allowable PA is the percentage of allowable aliasing energy.
Technically speaking, it the amount of aliasing which can be tolerated. In practice, a bigger PA will produce small localization error of the edges. I hope to give you a more practical explanation soon...

Then, you have to choose between three different LoG implementations:
* Standard LoG: The LoG convolution followed by a zero crossing of the resulting image
* LoG with Roberts: The Standard LoG with a Roberts gradient.
* LoG with Sobel: The Standard LoG with a Sobel gradient.

The gradient step is useful to skip spurious contours in the image. If you choose to apply a gradient, then you can control its thresholds. 

In practice, only data between the PC1% and the PC2% is kept after the gradient computation.

Again, the suggested values will work in most cases.

Finally, you must specify a standard deviation to apply to the LoG filter. In practice, the larger it is, the few details you'll get.

You can play with it but, in most cases, the suggested value of 2.0 will work quite well.

Now, a few examples...

This is a church in Seydisfjordur, Iceland.

![](/docs/logchurch.jpeg)

Now we apply a standard LoG filter with PA=0.00001 and standard deviation of 2.0
Look at how well the church's contours are outlined. But you also get a lot of edges in the background, caused by small differences in the gray levels.
Notice that all contours are closed, an interesting property of the LoG filter.

![](/docs/logchurch1.gif)

Undo this filter and apply a LoG filter with Sobel gradient convolution, same parameters as before, PC1=25, PC2=60.
The image is now cleaner than before, but we have also lost details about the church.

![](/docs/logchurch2.gif)

Undo again. Apply now a standard LoG filter with PA=0.00001 and standard deviation of 4.0
Now the church is roughly outlined.

![](/docs/logchurch3.gif)

I think you'll make an artistic use of this filter.
