/*
** This file is part of OSPREY 3.0
** 
** OSPREY Protein Redesign Software Version 3.0
** Copyright (C) 2001-2018 Bruce Donald Lab, Duke University
** 
** OSPREY is free software: you can redistribute it and/or modify
** it under the terms of the GNU General Public License version 2
** as published by the Free Software Foundation.
** 
** You should have received a copy of the GNU General Public License
** along with OSPREY.  If not, see <http://www.gnu.org/licenses/>.
** 
** OSPREY relies on grants for its development, and since visibility
** in the scientific literature is essential for our success, we
** ask that users of OSPREY cite our papers. See the CITING_OSPREY
** document in this distribution for more information.
** 
** Contact Info:
**    Bruce Donald
**    Duke University
**    Department of Computer Science
**    Levine Science Research Center (LSRC)
**    Durham
**    NC 27708-0129
**    USA
**    e-mail: www.cs.duke.edu/brd/
** 
** <signature of Bruce Donald>, Mar 1, 2018
** Bruce Donald, Professor of Computer Science
*/

package edu.duke.cs.osprey.tools;

import edu.duke.cs.osprey.structure.Atom;

/**
 *
 * @author mhall44
 */
public class PeriodicTable {
    //Gives information about chemical elements
    
    
    
    public static void setElementProperties(Atom atom, String elementType){
        //set atom to have the specified elementType, filling in its element-based properties accordingly
        setElementPropertiesByName(atom,elementType);
    }
    
    public static void inferElementProperties(Atom atom){
        //infer the element type of an atom from its name, and set its element-based properties accordingly
        setElementPropertiesByName(atom,atom.name);
    }
    
    
    
    
    //this function sets element-based properties (element type, atomic number,
    //radius, and mass) for an atom based on a name (optimally the name of the element,
    //but atom name if element name not available)
    private static void setElementPropertiesByName(Atom atom, String name){
        // Sets default properties such as radii, element number, mass, element type
        // RHL I don't know where the radii values come from but they are good estimates
        //  usually atoms will be assigned better information from a forcefield
    
            int elementNumber=0;
            double radius=0, mass=0;
            String elementType=null;

            if ( name.indexOf("Ac") != -1 ){
                    elementNumber = 89;
                    radius = 295;
                    mass = 227;
                    elementType = "Ac";
            }
            else if( name.indexOf("Ag") != -1 ){
                    radius = 398;
                    elementNumber = 47;
                    mass = 107.9;
                    elementType = "Ag";
            }
            else if( name.indexOf("Al") != -1 ){
                    radius = 338;
                    elementNumber = 13;
                    mass = 26.98;
                    elementType = "Al";
            }
            else if( name.indexOf("Am") != -1 ){
                    elementNumber = 95;
                    radius = 230;
                    mass = 243;
                    elementType = "Am";
            }
            else if( name.indexOf("Ar") != -1 ){
                    elementNumber = 18;
                    radius = 392;
                    mass = 39.95;
                    elementType = "Ar";
            }
            else if( name.indexOf("As") != -1 ){
                    elementNumber = 33;
                    radius = 302;
                    mass = 74.92;
                    elementType = "As";
            }
            else if( name.indexOf("At") != -1 ){
                    elementNumber = 85;
                    radius = 302;
                    mass = 210;
                    elementType = "At";
            }
            else if( name.indexOf("Au") != -1 ){
                    radius = 375;
                    elementNumber = 79;
                    mass = 197;
                    elementType = "Au";
            }
            else if( name.indexOf("Ba") != -1 ){
                    radius = 335;
                    elementNumber = 56;
                    mass = 137.3;
                    elementType = "Ba";
            }
            else if( name.indexOf("Be") != -1 ){
                    elementNumber = 4;
                    radius = 88;
                    mass = 9.012;
                    elementType = "Be";
            }
            else if( name.indexOf("Bi") != -1 ){
                    elementNumber = 83;
                    radius = 385;
                    mass = 209;
                    elementType = "Bi";
            }
            else if( name.indexOf("Bk") != -1 ){
                    elementNumber = 97;
                    radius = 225;
                    mass = 247;
                    elementType = "Bk";
            }
            else if( name.indexOf("Br") != -1 ){
                    radius = 302;
                    elementNumber = 97;
                    mass = 247;
                    elementType = "Br";
            }
            else if( name.indexOf("Ca") != -1 ){
                    radius = 248;
                    elementNumber = 20;
                    mass = 40.08;
                    elementType = "Ca";
            }
            else if( name.indexOf("Cd") != -1 ){
                    elementNumber = 48;
                    radius = 422;
                    mass = 112.4;
                    elementType = "Cd";
            }
            else if( name.indexOf("Ce") != -1 ){
                    elementNumber = 58;
                    radius = 458;
                    mass = 140.1;
                    elementType = "Ce";
            }
            else if( name.indexOf("Cf") != -1 ){
                    elementNumber = 98;
                    radius = 222;
                    mass = 249;
                    elementType = "Cf";
            }
            else if( name.indexOf("Cl") != -1 ){
                    elementNumber = 17;
                    radius = 250;
                    mass = 35.45;
                    elementType = "Cl";
            }
            else if( name.indexOf("Cm") != -1 ){
                    elementNumber = 96;
                    radius = 228;
                    mass = 247;
                    elementType = "Cm";
            }
            else if( name.indexOf("Co") != -1 ){
                    radius = 332;
                    elementNumber = 27;
                    mass = 58.93;
                    elementType = "Co";
            }
            else if( name.indexOf("Cr") != -1 ){
                    radius = 338;
                    elementNumber = 24;
                    mass = 52;
                    elementType = "Cr";
            }
            else if( name.indexOf("Cs") != -1 ){
                    elementNumber = 55;
                    radius = 418;
                    mass = 132.9;
                    elementType = "Cs";
            }
            else if( name.indexOf("Cu") != -1 ){
                    radius = 380;
                    elementNumber = 29;
                    mass = 63.55;
                    elementType = "Cu";
            }
            else if( name.indexOf("Dy") != -1 ){
                    elementNumber = 66;
                    radius = 438;
                    mass = 162.5;
                    elementType = "Dy";
            }
            else if( name.indexOf("Er") != -1 ){
                    elementNumber = 68;
                    radius = 432;
                    mass = 167.3;
                    elementType = "Er";
            }
            else if( name.indexOf("Es") != -1 ){
                    elementNumber = 99;
                    radius = 220;
                    mass = 254;
                    elementType = "Es";
            }
            else if( name.indexOf("Eu") != -1 ){
                    elementNumber = 63;
                    radius = 498;
                    mass = 152;
                    elementType = "Eu";
            }
            else if( name.indexOf("Fe") != -1 ){
                    radius = 335;
                    elementNumber = 26;
                    mass = 55.85;
                    elementType = "Fe";
            }
            else if( name.indexOf("Fm") != -1 ){
                    elementNumber = 100;
                    radius = 218;
                    mass = 250;
                    elementType = "Fm";
            }
            else if( name.indexOf("Fr") != -1 ){
                    elementNumber = 87;
                    radius = 450;
                    mass = 223;
                    elementType = "Fr";
            }
            else if( name.indexOf("Ga") != -1 ){
                    radius = 305;
                    elementNumber = 31;
                    mass = 69.72;
                    elementType = "Ga";
            }
            else if( name.indexOf("Gd") != -1 ){
                    elementNumber = 64;
                    radius = 448;
                    mass = 157.3;
                    elementType = "Gd";
            }
            else if( name.indexOf("Ge") != -1 ){
                    elementNumber = 32;
                    radius = 292;
                    mass = 72.59;
                    elementType = "Ge";
            }
            else if( name.indexOf("He") != -1 ){
                    radius = 400;
                    elementNumber = 2;
                    mass = 4.003;
                    elementType = "He";
            }
            else if( name.indexOf("Hf") != -1 ){
                    elementNumber = 72;
                    radius = 392;
                    mass = 178.5;
                    elementType = "Hf";
            }
            else if( name.indexOf("Hg") != -1 ){
                    elementNumber = 80;
                    radius = 425;
                    mass = 200.6;
                    elementType = "Hg";
            }
            else if( name.indexOf("Ho") != -1 ){
                    elementNumber = 67;
                    radius = 435;
                    mass = 164.9;
                    elementType = "Ho";
            }
            else if( name.indexOf("In") != -1 ){
                    elementNumber = 49;
                    radius = 408;
                    mass = 114.8;
                    elementType = "In";
            }
            else if( name.indexOf("Ir") != -1 ){
                    elementNumber = 77;
                    radius = 330;
                    mass = 192.2;
                    elementType = "Ir";
            }
            else if( name.indexOf("Kr") != -1 ){
                    elementNumber = 36;
                    radius = 400;
                    mass = 83.8;
                    elementType = "Kr";
            }
            else if( name.indexOf("La") != -1 ){
                    elementNumber = 57;
                    radius = 468;
                    mass = 138.9;
                    elementType = "La";
            }
            else if( name.indexOf("Li") != -1 ){
                    radius = 170;
                    elementNumber = 3;
                    mass = 6.941;
                    elementType = "Li";
            }
            else if( name.indexOf("Lr") != -1 ){
                    elementNumber = 103;
                    radius = 210;
                    mass = 257;
                    elementType = "Lr";
            }
            else if( name.indexOf("Lu") != -1 ){
                    elementNumber = 71;
                    radius = 430;
                    mass = 175;
                    elementType = "Lu";
            }
            else if( name.indexOf("Md") != -1 ){
                    elementNumber = 101;
                    radius = 215;
                    mass = 256;
                    elementType = "Md";
            }
            else if( name.indexOf("Mg") != -1 ){
                    radius = 275;
                    elementNumber = 12;
                    mass = 24.31;
                    elementType = "Mg";
            }
            else if( name.indexOf("Mn") != -1 ){
                    radius = 338;
                    elementNumber = 25;
                    mass = 54.94;
                    elementType = "Mn";
            }
            else if( name.indexOf("Mo") != -1 ){
                    elementNumber = 42;
                    radius = 368;
                    mass = 95.94;
                    elementType = "Mo";
            }
            else if( name.indexOf("Na") != -1 ){
                    elementNumber = 11;
                    radius = 243;
                    mass = 22.99;
                    elementType = "Na";
            }
            else if( name.indexOf("Nb") != -1 ){
                    elementNumber = 41;
                    radius = 370;
                    mass = 92.91;
                    elementType = "Nb";
            }
            else if( name.indexOf("Nd") != -1 ){
                    elementNumber = 60;
                    radius = 452;
                    mass = 144.2;
                    elementType = "Nd";
            }
            else if( name.indexOf("Ne") != -1 ){
                    elementNumber = 10;
                    radius = 280;
                    mass = 20.18;
                    elementType = "Ne";
            }
            else if( name.indexOf("Ni") != -1 ){
                    radius = 405;
                    elementNumber = 28;
                    mass = 58.69;
                    elementType = "Ni";
            }
            else if( name.indexOf("No") != -1 ){
                    elementNumber = 102;
                    radius = 212;
                    mass = 253;
                    elementType = "No";
            }
            else if( name.indexOf("Np") != -1 ){
                    elementNumber = 93;
                    radius = 238;
                    mass = 237;
                    elementType = "Np";
            }
            else if( name.indexOf("Os") != -1 ){
                    elementNumber = 76;
                    radius = 342;
                    mass = 190.2;
                    elementType = "Os";
            }
            else if( name.indexOf("Pa") != -1 ){
                    elementNumber = 91;
                    radius = 222;
                    mass = 231;
                    elementType = "Pa";
            }
            else if( name.indexOf("Pb") != -1 ){
                    elementNumber = 82;
                    radius = 385;
                    mass = 207.2;
                    elementType = "Pb";
            }
            else if( name.indexOf("Pd") != -1 ){
                    elementNumber = 46;
                    radius = 375;
                    mass = 106.4;
                    elementType = "Pd";
            }
            else if( name.indexOf("Pm") != -1 ){
                    elementNumber = 61;
                    radius = 450;
                    mass = 147;
                    elementType = "Pm";
            }
            else if( name.indexOf("Po") != -1 ){
                    elementNumber = 84;
                    radius = 420;
                    mass = 210;
                    elementType = "Po";
            }
            else if( name.indexOf("Pr") != -1 ){
                    elementNumber = 59;
                    radius = 455;
                    mass = 140.9;
                    elementType = "Pr";
            }
            else if( name.indexOf("Pt") != -1 ){
                    elementNumber = 78;
                    radius = 375;
                    mass = 195.1;
                    elementType = "Pt";
            }
            else if( name.indexOf("Pu") != -1 ){
                    elementNumber = 94;
                    radius = 232;
                    mass = 242;
                    elementType = "Pu";
            }
            else if( name.indexOf("Ra") != -1 ){
                    elementNumber = 88;
                    radius = 358;
                    mass = 226;
                    elementType = "Ra";
            }
            else if( name.indexOf("Rb") != -1 ){
                    elementNumber = 37;
                    radius = 368;
                    mass = 85.47;
                    elementType = "Rb";
            }
            else if( name.indexOf("Re") != -1 ){
                    elementNumber = 75;
                    radius = 338;
                    mass = 186.2;
                    elementType = "Re";
            }
            else if( name.indexOf("Rh") != -1 ){
                    elementNumber = 45;
                    radius = 362;
                    mass = 102.9;
                    elementType = "Rh";
            }
            else if( name.indexOf("Rn") != -1 ){
                    elementNumber = 88;
                    radius = 475;
                    mass = 222;
                    elementType = "Rn";
            }
            else if( name.indexOf("Ru") != -1 ){
                    elementNumber = 44;
                    radius = 350;
                    mass = 101.1;
                    elementType = "Ru";
            }
            else if( name.indexOf("Sb") != -1 ){
                    elementNumber = 51;
                    radius = 365;
                    mass = 121.8;
                    elementType = "Sb";
            }
            else if( name.indexOf("Sc") != -1 ){
                    radius = 360;
                    elementNumber = 21;
                    mass = 44.96;
                    elementType = "Sc";
            }
            else if( name.indexOf("Se") != -1 ){
                    elementNumber = 34;
                    radius = 305;
                    mass = 78.96;
                    elementType = "Se";
            }
            else if( name.indexOf("Si") != -1 ){
                    radius = 300;
                    elementNumber = 14;
                    mass = 28.09;
                    elementType = "Si";
            }
            else if( name.indexOf("Sm") != -1 ){
                    elementNumber = 62;
                    radius = 450;
                    mass = 150.4;
                    elementType = "Sm";
            }
            else if( name.indexOf("Sn") != -1 ){
                    elementNumber = 50;
                    radius = 365;
                    mass = 118.7;
                    elementType = "Sn";
            }
            else if( name.indexOf("Sr") != -1 ){
                    elementNumber = 38;
                    radius = 280;
                    mass = 87.62;
                    elementType = "Sr";
            }
            else if( name.indexOf("Ta") != -1 ){
                    elementNumber = 73;
                    radius = 358;
                    mass = 180.9;
                    elementType = "Ta";
            }
            else if( name.indexOf("Tb") != -1 ){
                    elementNumber = 65;
                    radius = 440;
                    mass = 158.9;
                    elementType = "Tb";
            }
            else if( name.indexOf("Tc") != -1 ){
                    elementNumber = 43;
                    radius = 338;
                    mass = 99;
                    elementType = "Tc";
            }
            else if( name.indexOf("Te") != -1 ){
                    elementNumber = 52;
                    radius = 368;
                    mass = 127.6;
                    elementType = "Te";
            }
            else if( name.indexOf("Th") != -1 ){
                    elementNumber = 90;
                    radius = 255;
                    mass = 232;
                    elementType = "Th";
            }
            else if( name.indexOf("Ti") != -1 ){
                    radius = 368;
                    elementNumber = 22;
                    mass = 47.88;
                    elementType = "Ti";
            }
            else if( name.indexOf("Tl") != -1 ){
                    elementNumber = 81;
                    radius = 388;
                    mass = 204.4;
                    elementType = "Tl";
            }
            else if( name.indexOf("Tm") != -1 ){
                    elementNumber = 69;
                    radius = 430;
                    mass = 168.9;
                    elementType = "Tm";
            }
            else if( name.indexOf("Une") != -1 ){
                    elementNumber = 109;
                    mass = 266;
                    elementType = "Une";
            }
            else if( name.indexOf("Unh") != -1 ){
                    elementNumber = 106;
                    mass = 263;
                    elementType = "Unh";
            }
            else if( name.indexOf("Uno") != -1 ){
                    elementNumber = 108;
                    mass = 265;
                    elementType = "Uno";
            }
            else if( name.indexOf("Unp") != -1 ){
                    elementNumber = 105;
                    mass = 260;
                    elementType = "Unp";
            }
            else if( name.indexOf("Unq") != -1 ){
                    elementNumber = 104;
                    mass = 257;
                    elementType = "Unq";
            }
            else if( name.indexOf("Uns") != -1 ){
                    elementNumber = 107;
                    mass = 262;
                    elementType = "Uns";
            }
            else if( name.indexOf("Xe") != -1 ){
                    elementNumber = 54;
                    radius = 425;
                    mass = 131.3;
                    elementType = "Xe";
            }
            else if( name.indexOf("Yb") != -1 ){
                    elementNumber = 70;
                    radius = 485;
                    mass = 173;
                    elementType = "Yb";
            }
            else if( name.indexOf("Zn") != -1 ){
                    radius = 362;
                    elementNumber = 30;
                    mass = 65.39;
                    elementType = "Zn";
            }
            else if( name.indexOf("Zr") != -1 ){
                    elementNumber = 40;
                    radius = 390;
                    mass = 91.22;
                    elementType = "Zr";
            }
            //The radii for C, H, N, O, P, S were changed to those used by the Richardsons' PROBE;
            //		The old radii are commented out to the side
            else if( name.indexOf("C") == 0 ){
                    radius = 165; //radius = 180;
                    elementNumber = 6;
                    mass = 12.01;
                    elementType = "C";
            }
            else if( (name.indexOf("H") == 0) || 
                            ( (name.indexOf("H") == 1) && ((name.charAt(0)>='0') && (name.charAt(0)<='9')) ) ){
                    radius = 100; //radius = 80;
                    elementNumber = 1;
                    mass = 1;
                    elementType = "H";
            }
            else if( name.indexOf("N") == 0 ){
                    radius = 155; //radius = 170;
                    elementNumber = 7;
                    mass = 14.01;
                    elementType = "N";
            }
            else if( name.indexOf("O") == 0 ){
                    radius = 140; //radius = 170;
                    elementNumber = 8;
                    mass = 16;
                    elementType = "O";
            }
            else if( name.indexOf("B") == 0 ){
                    radius = 208;
                    elementNumber = 5;
                    mass = 10.81;
                    elementType = "B";
            }
            else if( name.indexOf("I") == 0 ){
                    radius = 350;
                    elementNumber = 53;
                    mass = 126.9;
                    elementType = "I";
            }
            else if( name.indexOf("F") == 0 ){
                    radius = 160;
                    elementNumber = 9;
                    mass = 19.0;
                    elementType = "F";
            }
            else if( name.indexOf("P") == 0 ){
                    radius = 180; //radius = 259;
                    elementNumber = 15;
                    mass = 30.97;
                    elementType = "P";
            }
            else if( name.indexOf("K") == 0 ){
                    radius = 332;
                    elementNumber = 19;
                    mass = 39.1;
                    elementType = "K";
            }
            else if( name.indexOf("S") == 0 ){
                    radius = 180; //radius = 255;
                    elementNumber = 16;
                    mass = 32.07;
                    elementType = "S";
            }
            else if( name.indexOf("U") == 0 ){
                    radius = 242;
                    elementNumber = 92;
                    mass = 238;
                    elementType = "U";
            }
            else if( name.indexOf("V") == 0 ){
                    radius = 332;
                    elementNumber = 23;
                    mass = 50.94;
                    elementType = "V";
            }
            else if( name.indexOf("W") == 0 ){
                    radius = 342;
                    elementNumber = 74;
                    mass = 183.9;
                    elementType = "W";
            }
            else if( name.indexOf("Y") == 0 ){
                    radius = 445;
                    elementNumber = 39;
                    mass = 88.91;
                    elementType = "Y";
            }
            else {  //unrecognized atom type
                    radius = 10000;
                    elementType = "DU"; //dummy atom
            }
            
            //now actually set the properties in the atom
            atom.elementType = elementType;
            atom.radius = radius;
            atom.mass = mass;
            atom.elementNumber = elementNumber;
    }
    
    
}
