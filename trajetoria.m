function traj =  trajetoria(a,v_0, y_0, t_f)
    % tempo ate o primeiro acionamento
    tempo = [0:0.0001:t_f];

    % trajetoria em queda livre
    v1 = v_0 + a.*tempo;
    y1 = y_0 + v_0.*tempo + a.*tempo.^2/2;

    traj.v = v1;
    traj.y = y1;
end