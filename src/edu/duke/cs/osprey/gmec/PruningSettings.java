/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.gmec;

import java.io.Serializable;

/**
 *
 * @author mhall44
 */
public class PruningSettings implements Serializable {
    public boolean typedep = false;
    public double boundsThresh = 100;
    public int algOption = 1;
    public boolean useFlags = true;
    public boolean useTriples = false;
    public double stericThresh = 100;
}
