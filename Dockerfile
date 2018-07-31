FROM gcc:4.9

MAINTAINER Matthieu Foll <follm@iarc.fr>

WORKDIR /usr/src/myapp

RUN git clone https://github.com/mfoll/BayeScanHierachical.git && \
	mv BayeScanHierachical/* /usr/src/myapp && \
	make 

CMD ["./bayescanH"]
