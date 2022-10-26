function [ac_ModeS_ns, v_ModeS_ns, u_ModeS_nsi1] = ModeS(order, VV, diagM, ACC_eli, F_ModeS_ps, ps_weizhi, Ken, an, bn, qni)

    for ii = 1:order
        PACC_el = VV(:,ii)' * (ACC_eli * diagM + F_ModeS_ps * ps_weizhi);
    end
    %PPn = -PACC_el(:, i) - an * qn(: , i) - bn * qn(: , i-1);
    PPn = -PACC_el - an * qni - bn * qni_1;
    qni1=Ken \ PPn;
    
    for iiii = 1:order
        u_ModeS_nsi1 = u_ModeS_nsi1 + VV(:,iiii)'* qni1;
    end
    
    v_ModeS_ns = (u_ModeS_nsi1 - u_ModeS_nsi_1) / (dt*2);
    ac_ModeS_ns = (u_ModeS_nsi1 - 2 * u_ModeS_nsi + u_ModeS_nsi_1) / (dt^2);
end