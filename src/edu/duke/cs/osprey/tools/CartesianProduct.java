/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.tools;

/**
 *
 * @author hmn5
 */
import java.util.ArrayList;
import edu.duke.cs.osprey.tools.ObjectIO;
public class CartesianProduct<T> {
    //This provides general tools for taking the cartesian product of arbitary sets
    //implementd as arraylists of any type
    
    
    ArrayList<ArrayList<T>> originalList;
    ArrayList<ArrayList<T>> cartProduct;
    
    //We are given a set of sets implemented as arraylists
    //Given {A,B,C,...} we want {(a,b,c,..) | a in A, b in B, c in C...}
    public CartesianProduct (ArrayList<ArrayList<T>> originalSets){
        this.originalList = originalSets;
        if (originalSets.size()==2){
            //Given {A,B} we return {(a,b) | a in A and b in B}
            this.cartProduct = cartProductTwoSets(originalSets.get(0), originalSets.get(1));
        }
        else{
            //We are given {A,B,C,...}
            ArrayList<ArrayList<T>> productSet  = cartProductTwoSets(originalSets.get(0), originalSets.get(1));
            for (int setIndex = 2; setIndex<originalSets.size(); setIndex++){
                productSet = cartProductAddSet(productSet, originalSets.get(setIndex));
            }
            this.cartProduct = productSet;
        }
        }
    
    public ArrayList<ArrayList<T>> getCartProduct(){
        return this.cartProduct;
    }
    
    //returns cartesian product between two sets
    //{(a,b) | a in A and b in B}
    public static <T> ArrayList<ArrayList<T>> cartProductTwoSets(ArrayList<T> setA, ArrayList<T> setB){
        ArrayList<ArrayList<T>> productSet = new ArrayList<ArrayList<T>>();
        for (int aIndex = 0; aIndex<setA.size(); aIndex++){
            for (int bIndex = 0; bIndex<setB.size(); bIndex++){
                ArrayList<T> product = new ArrayList<T>();
                product.add(setA.get(aIndex));
                product.add(setB.get(bIndex));
                productSet.add(product);
            }
        }
        return productSet;
    }
    
    /*This is a utility function to compute cartesian products of multiple sets
     This will be an intermediate set where we have the cartesian product of sets A and B
     and we want to add one more set C. So we have {(a,b)|a in A and b in B} and we will add set C
     to get {(a,b,c) | a in A, b in B, c in C}
     */
    
    public static <T> ArrayList<ArrayList<T>> cartProductAddSet(ArrayList<ArrayList<T>> productSet, ArrayList<T> setC){
        ArrayList<ArrayList<T>> originalProductSet = productSet;
        ArrayList<ArrayList<T>> newProductSet = new ArrayList<ArrayList<T>>();
        for (T element : setC){
            for (ArrayList<T> productElement : originalProductSet){
                ArrayList<T> newElement = (ArrayList<T>) ObjectIO.deepCopy(productElement);
                newElement.add(element);
                newProductSet.add(newElement);
            }
        }
        return newProductSet;
    }
    
    //If we don't want to create a CartesianProduct Object, we can simply use this static method
    //It does the same thing as the constructor
    public static <T> ArrayList<ArrayList<T>> cartesianProduct(ArrayList<ArrayList<T>> originalSets){
        ArrayList<ArrayList<T>> productSet;
        if (originalSets.size() == 2){
            productSet = cartProductTwoSets(originalSets.get(0), originalSets.get(1));
        }
        else{
            productSet  = cartProductTwoSets(originalSets.get(0), originalSets.get(1));
            for (int setIndex = 2; setIndex<originalSets.size(); setIndex++){
                productSet = cartProductAddSet(productSet, originalSets.get(setIndex));
            }
        }
        return productSet;
    }
}
