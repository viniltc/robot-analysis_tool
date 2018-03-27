function runs = get_runs_dyad(data)

runs = data.state'==2 | data.state'==3 | data.state'==4;
%runs = data.state'>=2 & data.state'<=4;
