function runs = get_runs(data)

runs = data.state==2 | data.state==3; 
