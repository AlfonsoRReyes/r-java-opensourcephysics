/*
 * Open Source Physics software is free software as described near the bottom of this code file.
 *
 * For additional information and documentation on Open Source Physics please see: 
 * <http://www.opensourcephysics.org/>
 */

package org.opensourcephysics.manual.ch11;
import org.opensourcephysics.display3d.simple3d.*;
import org.opensourcephysics.frames.Display3DFrame;
import org.opensourcephysics.numerics.*;

/**
 * QRotationApp demonstrates how to use the Quaternion class to rotate an Element.
 *
 * @author W. Christian
 * @version 1.0
 */
public class QRotationApp {
  public static void main(String[] args) {
    Display3DFrame frame = new Display3DFrame("Quaternion Rotation");
    frame.setDecorationType(org.opensourcephysics.display3d.core.VisualizationHints.DECORATION_AXES);
    Element ellipsoid = new ElementEllipsoid();
    Element arrow = new ElementArrow();
    ellipsoid.setSizeXYZ(0.4, 0.4, 1.0);
    Quaternion rotation = new Quaternion(1, 0, 0, 0);
    ellipsoid.setTransformation(rotation);
    frame.addElement(arrow);
    frame.addElement(ellipsoid);
    frame.setVisible(true);
    frame.setDefaultCloseOperation(javax.swing.JFrame.EXIT_ON_CLOSE);
    double theta = 0;
    double[] n = new double[] {0.5, 0.5, 1}; // rotation axis
    arrow.setSizeXYZ(n);
    while(true) { // animate until the program exits
      try {
        Thread.sleep(100);
      } catch(InterruptedException ex) {}
      theta += Math.PI/40;
      double cos = Math.cos(theta/2), sin = Math.sin(theta/2);
      // setCoordinates will normalize quaternion
      rotation.setCoordinates(cos, sin*n[0], sin*n[1], sin*n[2]);
      ellipsoid.setTransformation(rotation);
      frame.render();
    }
  }
}

/*
 * Open Source Physics software is free software; you can redistribute
 * it and/or modify it under the terms of the GNU General Public License (GPL) as
 * published by the Free Software Foundation; either version 2 of the License,
 * or(at your option) any later version.

 * Code that uses any portion of the code in the org.opensourcephysics package
 * or any subpackage (subdirectory) of this package must must also be be released
 * under the GNU GPL license.
 *
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston MA 02111-1307 USA
 * or view the license online at http://www.gnu.org/copyleft/gpl.html
 *
 * For additional information and documentation on Open Source Physics,
 * please see <http://www.opensourcephysics.org/>.
 *
 * Copyright (c) 2007  The Open Source Physics project
 *                     http://www.opensourcephysics.org
 */
