# Django-Chemical_Engineering
Django project to perform some Process Engineering concepts on a website.

Originally developed under Django 3.0.5 environment with Python 3.8

The VLE application is used to perform Vapor-Liquid Equilibrium calculations. 
For now only Fugacity and Activity coefficient calculations are complete to serve directly on a website. 

After you have created a Django Project for yourself, run the "get_chemical()" command on VLE/chemsep_operation file.
This will parse the XML file and then create Python class pickle objects in VLE/chemicals, which are used as chemical compounds. As long as the file structure and "class Chemical()" definition in the "VLE/chemsep_operation" file remain unchanged, pickle objects are created only once.

VLE part of this project is live on the website www.chemicaleasy.com 
