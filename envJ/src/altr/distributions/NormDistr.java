/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package altr.distributions;

import java.util.Random;

/**
 *
 * @author Vitaly
 */
public class NormDistr implements Distribution {
    
    private final Random generator = new Random();
    private final double mean;
    private final double deviation;
    
    
    public NormDistr(double mean, double deviation) {
        this.mean = mean;
        this.deviation = deviation;
    }

    @Override
    public double getValue() {
        return this.mean + (this.deviation * generator.nextGaussian());
    }

    @Override
    public String getName() {
        return "NormDistr";
    }
    
}
