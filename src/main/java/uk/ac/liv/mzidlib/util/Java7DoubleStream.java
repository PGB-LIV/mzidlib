/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package uk.ac.liv.mzidlib.util;

import java.util.Collection;

/**
 *
 * @author SPerkins
 */
public class Java7DoubleStream {
    private Collection<Double> collection = null;
    
    
    public static Java7DoubleStream doubleStream(Collection<Double> collection) {
        return new Java7DoubleStream(collection);
    }
    
    private Java7DoubleStream(Collection<Double> collection) {
        this.collection = collection;
    }
}
