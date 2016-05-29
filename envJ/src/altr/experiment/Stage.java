/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package altr.experiment;

import altr.distributions.Distribution;

/**
 *
 * @author Vitaly
 */
public class Stage {
    private long experimentId;
    private long stepNum;
    private Distribution distribution;
    private long order;

    public long getExperimentId() {
        return experimentId;
    }

    public void setExperimentId(long experimentId) {
        this.experimentId = experimentId;
    }

    public long getStepNum() {
        return stepNum;
    }

    public void setStepNum(long stepNum) {
        this.stepNum = stepNum;
    }

    public Distribution getDestribution() {
        return distribution;
    }

    public void setDestribution(Distribution destribution) {
        this.distribution = destribution;
    }

    public long getOrder() {
        return order;
    }

    public void setOrder(long order) {
        this.order = order;
    }

    public Stage(long experimentId, long stepNum, Distribution distribution, long order) {
        this.experimentId = experimentId;
        this.stepNum = stepNum;
        this.distribution = distribution;
        this.order = order;
    }   

    @Override
    public String toString() {
        return "Stage{" + "experimentId=" + experimentId + ", stepNum=" + stepNum + ", distribution=" + distribution.getName() + ", order=" + order + '}';
    }
    
}
