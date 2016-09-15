package system;

public class DataTransmission {
	
	private int type = 0; // 0: source data transmission. 1: intermediate data transmission. 
	
	private Sample sample;
	
	private Double transmittedVolume = 0d;
	
	private double bandwidth;
		
	public DataTransmission(int type, Sample sample, double transmittedVolume, double bandwidth){
		this.setType(type);
		this.setSample(sample);
		this.setBandwidth(bandwidth);
		this.setTransmittedVolume(transmittedVolume);
	}

	public void setType(int type) {
		this.type = type;
	}

	public int getType() {
		return type;
	}


	public void setBandwidth(double bandwidth) {
		this.bandwidth = bandwidth;
	}

	public double getBandwidth() {
		return bandwidth;
	}

	public Double getTransmittedVolume() {
		return transmittedVolume;
	}

	public void setTransmittedVolume(Double transmittedVolume) {
		this.transmittedVolume = transmittedVolume;
	}

	public Sample getSample() {
		return sample;
	}

	public void setSample(Sample sample) {
		this.sample = sample;
	}
}
