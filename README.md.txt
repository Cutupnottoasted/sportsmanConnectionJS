Cody Hoang
5/25/2023
codynguyen.hoang@gmail.com

Thanks to Hiroyuki Sato, Vladimmir Agafonkin and Tim A. Pastva!
"If I have seen further, it is by standing on the should of giants." - Sir Issac Newton

Javascript scripts used at Sportsman's Connection.

simplify.js - a translation of the original javascript by Vladimir Agafonkin
    to function inside adobe illustrator. (extracted from Leaflet) 
    This script is used to reduce the number of anchors in each pathitem.
    It is able to reduce the number of points better than the built in simplify.
    However, it only out performs the simplify tool when the pathitem's anchor count is over ~500
    The script is CPU intensive and iterates through the pathitem multiple times.
    the current goal is to improve it's performance time. 

bezierFit.js - a translation of 1998 Pastva, Tim A. study on bezier curve fitting
    By implementing Berstein's Polynomials and Guass-Newton's method, we are able to approximate
    a best fit bezier curve give a set of data points. By doing so, we are able to follow a set of 
    data points and fit a quadratic bezier curve to the given points. 

By implementing both simplify.js and bezierFit.js into a single script you can optimize the reduction
    of anchors in a pathitem and maintain the original integrity of the said pathitem.

Notes: simplify.js is currently only written to use the high_quality boolean flag 
    and bezierFit.js still requires mxbern2.m and grad7.m. (aff_angle.m nearly finished)

Resources:
https://github.com/shspage 
https://github.com/mourner/simplify-js
https://upload.wikimedia.org/wikipedia/commons/1/15/Bezier_curve_fitting_%28IA_beziercurvefitti1094532767%29.pdf

