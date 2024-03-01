from numpy import isnan, nan
def task(catch_i,src,flows_df):
    catch = src['CatchId'].unique()[catch_i]
    catch_src = src[src['CatchId']==catch].reset_index()
    flow = flows_df[flows_df['feature_id']==catch]['streamflow']
    if len(flow) != 1:
        flow = 0.0
    flow = float(flow)
    if flow == 0.0:
        flow_stage = 0.0
        catch_src.loc[0,'hand_flow'] = flow
        catch_src.loc[0,'hand_stage'] = flow_stage
        catch_src = catch_src[catch_src['hand_stage'] != '']
        return(catch_src)
    
    if isnan(flow):
        flow_stage = nan
        catch_src.loc[0,'hand_flow'] = flow
        catch_src.loc[0,'hand_stage'] = flow_stage
        catch_src = catch_src[catch_src['hand_stage'] != '']
        return(catch_src)

    else:
        for src_flow_i in range(len(catch_src['Discharge (m3s-1)'])):
            max_i = max(range(len(catch_src['Discharge (m3s-1)'])))
            src_flow = catch_src['Discharge (m3s-1)'][src_flow_i]
            if max(catch_src['Discharge (m3s-1)']) == 0.0:
                catch_src.loc[src_flow_i,'hand_flow'] = flow
                catch_src.loc[src_flow_i,'hand_stage'] = 0.0
                catch_src = catch_src[catch_src['hand_stage'] != '']
                return(catch_src)
            if isnan(src_flow):
                src_flow = 0.0
                catch_src.loc[src_flow_i,'Discharge (m3s-1)'] = src_flow
            if src_flow_i == max_i:
                catch_src.loc[src_flow_i,'hand_flow'] = flow
                catch_src.loc[src_flow_i,'hand_stage'] = catch_src['Stage'][src_flow_i]
                catch_src = catch_src[catch_src['hand_stage'] != '']
                return(catch_src)
            if src_flow <= flow:
                continue
            else:
                flow_stage = catch_src['Stage'][src_flow_i-1] + (flow - catch_src['Discharge (m3s-1)'][src_flow_i - 1])*(catch_src['Stage'][src_flow_i] - catch_src['Stage'][src_flow_i-1])/(catch_src['Discharge (m3s-1)'][src_flow_i] - catch_src['Discharge (m3s-1)'][src_flow_i-1])
                catch_src.loc[src_flow_i,'hand_flow'] = flow
                catch_src.loc[src_flow_i,'hand_stage'] = flow_stage
                break
        catch_src = catch_src[catch_src['hand_stage'] != '']
        return(catch_src)
    
    