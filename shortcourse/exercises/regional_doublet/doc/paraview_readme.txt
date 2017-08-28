#Description: Instructions for using Paraview to create an image and a movie from stochastic_regional_doublet.h5.

0. Download Paraview executable from paraview.org. (It's free.)

1. Open stochastic_regional_doublet.h5 in Paraview
In Paraview, click on the open file icon.
Navigate file trees in the pop-up that opens to choose stochastic_regional_doublet.h5. 
An "Open Data With..." pop-up opens. Choose "PFLOTRAN Files"

2. Select data to load.
In the Properties pane, check the "Cell Arrays" box to select all output data. The "mesh" box is already checked. 
You don't need "Materials". 
Click "Apply". An outline of the model domain appears.

3. Setup model domain for viewing using the toolbar.
Click the +Y button.
From the dropdown menu that says "Outline" select "Surface".
From the dropdown menu that says "Solid Color" select "Kludged_material_ids_for_Visit"

4. Add vertical exaggeration using the Filters menu.
From the Filters menu choose "Alphabetical" and then "Transform".
In the Properties pane unclick "Show Box".
In the third box to right of "Scale" enter 5 (for 5x vertical exaggeration).
Click "Apply".
The transformed domain appears, and the original domain is removed from view. You may need to choose "Surface"
and/or "Kludged_material_ids_for_Visit" for the transformed domain.

5. Confirm you are looking at the transformed domain.
In the Pipeline Browser Pane, the eye icon is dark next to "Transform1" and grey next to "stochastic_regional_doublet.h5",
indicating that "Transform1" is visible and the other is not.

6. Split the screen.
Next to "RenderView1" in the upper right corner, click on the box that is divided with a horizontal line.
The viewing pane splits into top and bottom halves. In the new (bottom) half, choose "Render View".

7. Draw "Transform 1" in the bottom viewing pane.
The bottom (empty) viewing pane should be selected, indicated by it's blue outline. If it isn't selected, click
in it to select it.
In the Pipeline Browser pane, both eyes are grey. Click on the eye next to "Transform1" to display it.
Choose "Permeability_X [m^2]" from the left dropdown menu on the toolbar.
If necessary, choose "Surface" from the right dropdown menu on the toolbar.

8. Add camera link. (Assumes the bottom pane is selected with blue outline.)
From the Tools menu, choose "Add Camera Link..."
To create the link, follow the instructions in the pop-up that appears and click on the top pane.

9. Adjust the position of the model domain.
Use the left mouse button to rotate the domain.
Use the right mouse button to zoom.
Use the middle mouse button (or Shift + right mouse button) to pan.
Because they are linked, both views adjust simultaneously.

10. Remove the top two materials from view so you can see the heterogeneous permeability in Material 2.
With bottom viewing pane selected:
In the Pipeline Browser pane, click on "Transform1" to select it.
From the Filters menu choose "Common" and "Threshold".
In the Properties pane, choose "Kludged_material_ids_for_Visit" from the dropdown menu near the top.
Enter 1 in the Minimum box and 2 in the Maximum box.
Click "Apply".

11. Color the bottom (threshold) view by Tracer2 concentration.
With the bottom viewing pane and "Threshold1" selected, choose "Total_Tracer2 [M]" from the left dropdown
menu on the toolbar.

12. Color the top (transform) view by Tracer concentration.
Click at the edge of the top viewing pane to select it.
Click on "Transform1" in the Pipeline Browser to select it.
Choose "Total_Tracer [M] from the left dropdown menu on the toolbar.

13. Press the play button at the top of the screen to see concentration change with time.

14. Adjust the color scale in each viewing pane.
With top pane and "Transform1" selected, click on the double arrow icon at the left of the toolbar.
It autoscales to the current concentration range.
In the Color Map Editor pane, select "Use log scale when mapping data to colors".
Select the bottom pane and "Threshold1" and repeat.

15. Annotate time.
With the bottom pane selected:
From the Filters menu select "Annotation" and "Annotate Time Filter".
In the Properties pane, find "Time: %f" and edit to "Time: %0.2f y".
Click "Apply".

16. Add a title.
Select the top pane.
From the Sources menu choose "Text".
In the Properties pane, find "Text" and replace with "stochastic_regional_doublet".
Click "Apply".
(Do the same over again in the bottom pane to add "5x vertical exaggeration".)

17. Drag text, time, and color scales to desired positions.

18. Save all the work you've done setting this thing up.
From the File menu choose "Save State..."
In the pop-up that appears, navigate to desired location.
Type a file name ("st_reg_doublet") without extension.
Click "OK"

19. Save a movie
From the File menu choose "Save Animation..."
In the pop-up change 15.00 to 2.00 in the "Frame Rate" box.
Click "Save Animation"
In the next pop-up navigate to desired location, type a file name, and choose "AVI files (*.avi)".
Click "OK".

20. Save a screen shot.
At any time in the above sequence you can save a screen shot.
From the File menu choose "Save Screenshot..."
In the pop-up select "Save only selected view" if you want only the selected view.
(Otherwise you get both view panes in one image.)
Click "OK"
In the next pop-up navigate to desired location, type a file name, and click "OK".



