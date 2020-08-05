% populate binary matrix (m_allow) with 1 if ROI pair produce significant
% difference to random model in either comparison listed in (cmprsnnames).
function [m_allow] = Reduce_By_Rand_Signif(cmprsnnames,ad,mci,hc)

m_allow=zeros(132,132);
for kk = 1 : 2
    if cmprsnnames{kk}=='AD'
        m_allow=m_allow+ad;
    elseif cmprsnnames{kk}=='MC'
        m_allow=m_allow+mci;
    elseif cmprsnnames{kk}=='HC'
        m_allow=m_allow+hc;
    end
end
end