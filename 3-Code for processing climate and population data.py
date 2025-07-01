#人口数据改为用landscan，计算每个城市不同cdd范围内的人口数量  次国家尺度 2010 2019  2021 人口和cdd 都是1km  arcgis pro
#2010年和2021年使用不同的城市边界
# -*- coding: utf-8 -*-
import sys
import arcpy
import os
import pandas as pd
import numpy as np

# 设置输出编码为UTF-8，确保可以正确显示中文字符
sys.stdout.reconfigure(encoding='utf-8')

# 设置环境变量
arcpy.env.overwriteOutput = True
arcpy.env.scratchWorkspace = "L:/temp_workspace10"  # 确保路径存在

# 输入文件路径
cdd_raster_path = "M:/2019maxtem/2015-2020data/jisuan18.3_0.1_pop_year_0.125_SUM_change/0ssp126_rasters_1km/re_ssp126_2021.tif"
pop_raster_path = "F:/人口数据/lanscan_origin/landscan-global-2021-assets/landscan-global-2021_360.tif"
vector_path = "M:/arcgis_figs3_cityboundary/GHS_UCDB_MTUC_GLOBE_R2024A_GHSL_UCDB_MTUC_2020_GLOBE_R2024_2.shp"

# 输出文件路径
output_excel = "M:/2019maxtem/2015-2020data/jisuan18.3_0.1_pop_year_0.125_SUM_change/0ssp126_rasters_1km/CDD_POP-cities-21pop-21clim-21urban11.xlsx"

# 定义CDD区间，使用一个很大的数代替无穷大
cdd_bins = [0, 600, 1200, 1800, 2400, 3000, 1e10]
cdd_labels = ["lt600", "600_1200", "1200_1800", "1800_2400", "2400_3000", "gte3000"]

# 创建临时工作空间
temp_workspace = "L:/temp_workspace10"
if not os.path.exists(temp_workspace):
    os.makedirs(temp_workspace)
arcpy.env.workspace = temp_workspace

# 重分类CDD栅格到指定区间
reclass_ranges = arcpy.sa.RemapRange([[cdd_bins[i], cdd_bins[i + 1], i + 1] for i in range(len(cdd_bins) - 1)])
try:
    reclass_cdd = arcpy.sa.Reclassify(cdd_raster_path, "VALUE", reclass_ranges, "NODATA")
    reclass_cdd.save(os.path.join(temp_workspace, "reclass_cdd.tif"))
    print("重分类成功")
except arcpy.ExecuteError as e:
    print("重分类时出错: {}".format(e))
    raise

# 以每个国家为单元，对重分类后的CDD栅格和人口栅格进行区域统计
result_df = pd.DataFrame(columns=["Country"] + cdd_labels)

for cdd_label in cdd_labels:
    cdd_zone = cdd_labels.index(cdd_label) + 1
    zone_pop_raster = arcpy.sa.Con(reclass_cdd == cdd_zone, pop_raster_path)
    zonal_table = os.path.join(temp_workspace, "zonal_table_{}.dbf".format(cdd_label))
    try:
        arcpy.sa.ZonalStatisticsAsTable(vector_path, "ID_UC_G0_1", zone_pop_raster, zonal_table, "DATA", "SUM")
    except arcpy.ExecuteError as e:
        print(f"区域统计时出错 (CDD范围: {cdd_label}): {e}")
        continue

    # 读取统计结果并更新DataFrame
    stats_data = arcpy.da.TableToNumPyArray(zonal_table, ["ID_UC_G0_1", "SUM"])
    for row in stats_data:
        country_name = row['ID_UC_G0_1'].strip()
        pop_sum = row['SUM'] if row['SUM'] is not None else 0

        if country_name in result_df["Country"].values:
            result_df.loc[result_df["Country"] == country_name, cdd_label] = pop_sum
        else:
            new_row = {"Country": country_name, cdd_label: pop_sum}
            result_df = result_df.append(new_row, ignore_index=True)

# 保存结果到Excel
try:
    result_df.to_excel(output_excel, index=False)
    print("统计完成，结果保存在: {}".format(output_excel))
except Exception as e:
    print(f"保存结果到Excel时出错: {e}")

# 清理临时文件
try:
    arcpy.Delete_management(temp_workspace)
    print("临时文件已清理")
except arcpy.ExecuteError as e:
    print(f"清理临时文件时出错: {e}")





#2010年和2021年使用不同的城市边界，使用更大的城市边界
# -*- coding: utf-8 -*-
import sys
import arcpy
import os
import pandas as pd
import numpy as np

# 设置输出编码为UTF-8
sys.stdout.reconfigure(encoding='utf-8')

# 设置 ArcPy 环境
arcpy.env.overwriteOutput = True
arcpy.env.scratchWorkspace = "L:/temp_workspace10"

# 输入路径
cdd_raster_path = "M:/2019maxtem/2015-2020data/jisuan18.3_0.1_pop_year_0.125_SUM_change/0ssp126_rasters_1km/re_ssp126_2010.tif"
pop_raster_path = "F:/人口数据/lanscan_origin/landscan-global-2010-assets/landscan-global-2010_360.tif"
vector_path = "M:/arcgis_figs3_cityboundary/ucdb_urbanareas1.shp"

# 输出路径
output_excel = "M:/2019maxtem/2015-2020data/jisuan18.3_0.1_pop_year_0.125_SUM_change/0ssp126_rasters_1km/CDD_POP-cities-10pop-10clim-10urban-max1.xlsx"

# 设置CDD区间
cdd_bins = [0, 600, 1200, 1800, 2400, 3000, 1e10]
cdd_labels = ["lt600", "600_1200", "1200_1800", "1800_2400", "2400_3000", "gte3000"]

# 创建临时工作目录
temp_workspace = "L:/temp_workspace10"
if not os.path.exists(temp_workspace):
    os.makedirs(temp_workspace)
arcpy.env.workspace = temp_workspace

# 重分类CDD
reclass_ranges = arcpy.sa.RemapRange([[cdd_bins[i], cdd_bins[i + 1], i + 1] for i in range(len(cdd_bins) - 1)])
try:
    reclass_cdd = arcpy.sa.Reclassify(cdd_raster_path, "VALUE", reclass_ranges, "NODATA")
    reclass_cdd.save(os.path.join(temp_workspace, "reclass_cdd.tif"))
    print("重分类成功")
except arcpy.ExecuteError as e:
    print("重分类时出错: {}".format(e))
    raise

# 初始化结果数据
rows = []

# 每个区间分别进行Zonal统计
for cdd_label in cdd_labels:
    cdd_zone = cdd_labels.index(cdd_label) + 1
    zone_pop_raster = arcpy.sa.Con(reclass_cdd == cdd_zone, pop_raster_path)
    zonal_table = os.path.join(temp_workspace, f"zonal_table_{cdd_label}.dbf")

    try:
        arcpy.sa.ZonalStatisticsAsTable(vector_path, "name", zone_pop_raster, zonal_table, "DATA", "SUM")
    except arcpy.ExecuteError as e:
        print(f"区域统计出错（{cdd_label}）：{e}")
        continue

    # 读取统计结果并整理成字典形式
    stats_data = arcpy.da.TableToNumPyArray(zonal_table, ["name", "SUM"])
    for row in stats_data:
        country = row["name"].strip()
        pop_sum = row["SUM"] if row["SUM"] is not None else 0

        # 判断是否已有该城市的记录
        existing = next((r for r in rows if r["Country"] == country), None)
        if existing:
            existing[cdd_label] = pop_sum
        else:
            new_entry = {"Country": country, cdd_label: pop_sum}
            rows.append(new_entry)

# 转为DataFrame
result_df = pd.DataFrame(rows)
result_df = result_df.fillna(0)
result_df = result_df[["Country"] + cdd_labels]  # 确保列顺序

# 保存结果
try:
    result_df.to_excel(output_excel, index=False)
    print(f"统计完成，结果保存在: {output_excel}")
except Exception as e:
    print(f"保存结果时出错: {e}")

# 清理临时文件夹（可选）
try:
    arcpy.Delete_management(temp_workspace)
    print("临时文件已清理")
except arcpy.ExecuteError as e:
    print(f"清理临时文件时出错: {e}")






#新的代码
#人口数据改为用landscan，计算每个城市不同cdd范围内的人口数量  次国家尺度 2010 2019  2021 人口和cdd 都是1km  arcgis pro
#2010年和2021年使用不同的城市边界
# -*- coding: utf-8 -*-
import sys
import arcpy
import os
import pandas as pd
import numpy as np

# 设置输出编码为UTF-8，确保可以正确显示中文字符
sys.stdout.reconfigure(encoding='utf-8')

# 设置环境变量
arcpy.env.overwriteOutput = True
arcpy.env.scratchWorkspace = "L:/temp_workspace10"  # 确保路径存在

# 输入文件路径
cdd_raster_path = "M:/2019maxtem/2015-2020data/jisuan18.3_0.1_pop_year_0.125_SUM_change/0ssp126_rasters_1km/re_ssp126_2021.tif"
pop_raster_path = "F:/人口数据/lanscan_origin/landscan-global-2010-assets/landscan-global-2010_360.tif"
vector_path = "M:/arcgis_figs3_cityboundary/GHS_UCDB_MTUC_GLOBE_R2024A_GHSL_UCDB_MTUC_2010_GLOBE_R2024_2.shp"

# 输出文件路径
output_excel = "M:/2019maxtem/2015-2020data/jisuan18.3_0.1_pop_year_0.125_SUM_change/0ssp126_rasters_1km/CDD_POP-cities-10pop-21clim-10urban11.xlsx"

# 定义CDD区间，使用一个很大的数代替无穷大
cdd_bins = [0, 600, 1200, 1800, 2400, 3000, 1e10]
cdd_labels = ["lt600", "600_1200", "1200_1800", "1800_2400", "2400_3000", "gte3000"]

# 创建临时工作空间
temp_workspace = "L:/temp_workspace10"
if not os.path.exists(temp_workspace):
    os.makedirs(temp_workspace)
arcpy.env.workspace = temp_workspace

# 重分类CDD栅格到指定区间
reclass_ranges = arcpy.sa.RemapRange([[cdd_bins[i], cdd_bins[i + 1], i + 1] for i in range(len(cdd_bins) - 1)])
try:
    reclass_cdd = arcpy.sa.Reclassify(cdd_raster_path, "VALUE", reclass_ranges, "NODATA")
    reclass_cdd.save(os.path.join(temp_workspace, "reclass_cdd.tif"))
    print("重分类成功")
except arcpy.ExecuteError as e:
    print("重分类时出错: {}".format(e))
    raise

# 以每个国家为单元，对重分类后的CDD栅格和人口栅格进行区域统计
result_df = pd.DataFrame(columns=["Country"] + cdd_labels)

for cdd_label in cdd_labels:
    cdd_zone = cdd_labels.index(cdd_label) + 1
    zone_pop_raster = arcpy.sa.Con(reclass_cdd == cdd_zone, pop_raster_path)
    zonal_table = os.path.join(temp_workspace, "zonal_table_{}.dbf".format(cdd_label))
    try:
        arcpy.sa.ZonalStatisticsAsTable(vector_path, "ID_UC_G0_1", zone_pop_raster, zonal_table, "DATA", "SUM")
    except arcpy.ExecuteError as e:
        print(f"区域统计时出错 (CDD范围: {cdd_label}): {e}")
        continue

    # 读取统计结果并更新DataFrame
    stats_data = arcpy.da.TableToNumPyArray(zonal_table, ["ID_UC_G0_1", "SUM"])
    for row in stats_data:
        country_name = row['ID_UC_G0_1'].strip()
        pop_sum = row['SUM'] if row['SUM'] is not None else 0

        if country_name in result_df["Country"].values:
            result_df.loc[result_df["Country"] == country_name, cdd_label] = pop_sum
        else:
            new_row = pd.DataFrame([{"Country": country_name, cdd_label: pop_sum}])
            result_df = pd.concat([result_df, new_row], ignore_index=True)

# 保存结果到Excel
try:
    result_df.to_excel(output_excel, index=False)
    print("统计完成，结果保存在: {}".format(output_excel))
except Exception as e:
    print(f"保存结果到Excel时出错: {e}")

# 清理临时文件
try:
    arcpy.Delete_management(temp_workspace)
    print("临时文件已清理")
except arcpy.ExecuteError as e:
    print(f"清理临时文件时出错: {e}")









#统计次国家边界下城市面上的人口数据
import geopandas as gpd
import rasterio
from rasterstats import zonal_stats
import pandas as pd
import os

# ===== 路径配置 =====
subnational_path = r"M:\次国家行政边界数据\次国家边界1.shp"
city_path = r"M:\arcgis_figs3_cityboundary\GHS_UCDB_MTUC_GLOBE_R2024A_GHSL_UCDB_MTUC_2010_GLOBE_R2024_2.shp"
population_tif = r"F:\人口数据\lanscan_origin\landscan-global-2010-assets\landscan-global-2010_360.tif"
output_excel = r"M:\arcgis_figs3_cityboundary\次国家城市人口数据2010.xlsx"

# ===== 数据加载 =====
print("加载矢量数据...")
gdf_sub = gpd.read_file(subnational_path)
gdf_city = gpd.read_file(city_path)

# 坐标系检查与统一
if gdf_sub.crs != gdf_city.crs:
    print("城市边界坐标系与次国家不一致，正在转换...")
    gdf_city = gdf_city.to_crs(gdf_sub.crs)

# 空间索引加速
city_sindex = gdf_city.sindex

# ===== 核心处理逻辑 =====
print("开始统计每个次国家区域内城市人口数据...")

results = []

with rasterio.open(population_tif) as src:
    for idx, sub_row in gdf_sub.iterrows():
        fid_2 = sub_row['fid_2']  # 你关心的字段
        sub_geom = sub_row.geometry

        # 快速筛选与当前次国家区域相交的城市
        possible_index = list(city_sindex.intersection(sub_geom.bounds))
        possible_cities = gdf_city.iloc[possible_index]
        intersected_cities = possible_cities[possible_cities.intersects(sub_geom)]

        sub_pop_sum = 0  # 初始化当前次国家区域人口

        for _, city_row in intersected_cities.iterrows():
            inter_geom = city_row.geometry.intersection(sub_geom)
            if inter_geom.is_empty:
                continue

            stats = zonal_stats(
                inter_geom,
                population_tif,
                stats=["sum"],
                nodata=0,
                geojson_out=False
            )

            pop_sum = stats[0]["sum"] if stats and stats[0]["sum"] is not None else 0
            sub_pop_sum += pop_sum

        results.append({
            "fid_2": fid_2,
            "pop_sum_2010": sub_pop_sum
        })

# ===== 输出结果 =====
print("正在输出到 Excel...")
df_result = pd.DataFrame(results)
df_result.to_excel(output_excel, index=False)
print(f"完成！结果保存为：{output_excel}")




import geopandas as gpd
import rasterio
from rasterstats import zonal_stats
import pandas as pd
import warnings

warnings.filterwarnings("ignore")  # 忽略简化/拓扑等警告

# ========== 输入路径 ==========
subnat_shp_path = r"M:\次国家行政边界数据\次国家边界1.shp"
city_shp_path = r"M:\arcgis_figs3_cityboundary\GHS_UCDB_MTUC_GLOBE_R2024A_GHSL_UCDB_MTUC_2010_GLOBE_R2024_2.shp"
pop_raster_path = r"F:\人口数据\lanscan_origin\landscan-global-2010-assets\landscan-global-2010_360.tif"

# ========== 输出路径 ==========
output_excel_path = r"M:\arcgis_figs3_cityboundary\次国家城市人口数据2010_1.xlsx"

# ========== 步骤1：读取矢量数据 ==========
print("正在加载矢量数据...")
subnat_gdf = gpd.read_file(subnat_shp_path)
city_gdf = gpd.read_file(city_shp_path)

# ========== 步骤2：统一坐标系（强制为WGS84） ==========
target_crs = "EPSG:4326"
subnat_gdf = subnat_gdf.to_crs(target_crs)
city_gdf = city_gdf.to_crs(target_crs)

# ========== 步骤3：移除非法/空几何 ==========
city_gdf = city_gdf[city_gdf.is_valid & city_gdf.geometry.notnull()]
subnat_gdf = subnat_gdf[subnat_gdf.is_valid & subnat_gdf.geometry.notnull()]

# ========== 步骤4：遍历次国家区域统计城市人口 ==========
results = []

print("正在处理每个次国家区域...")
for idx, subnat_row in subnat_gdf.iterrows():
    fid_2 = subnat_row['fid_2']
    subnat_geom = subnat_row.geometry

    # 选取相交城市
    candidate_cities = city_gdf[city_gdf.geometry.intersects(subnat_geom)]

    if candidate_cities.empty:
        results.append({"fid_2": fid_2, "pop_sum": 0})
        continue

    # 裁剪城市为当前次国家内部部分
    subnat_geom_gdf = gpd.GeoDataFrame(geometry=[subnat_geom], crs=subnat_gdf.crs)
    try:
        clipped_cities = gpd.overlay(candidate_cities, subnat_geom_gdf, how='intersection')
    except Exception as e:
        print(f"fid_2={fid_2} overlay失败：{e}")
        results.append({"fid_2": fid_2, "pop_sum": 0})
        continue

    if clipped_cities.empty:
        results.append({"fid_2": fid_2, "pop_sum": 0})
        continue

    # 移除非法几何（保险）
    clipped_cities = clipped_cities[clipped_cities.is_valid & clipped_cities.geometry.notnull()]

    if clipped_cities.empty:
        results.append({"fid_2": fid_2, "pop_sum": 0})
        continue

    # 人口统计（避免nodata污染）
    try:
        stats = zonal_stats(
            clipped_cities,
            pop_raster_path,
            stats="sum",
            nodata=-9999,
            masked=True
        )
        # 防止负值污染：所有负值变为0
        total_pop = sum(max(s.get("sum", 0) or 0, 0) for s in stats)
    except Exception as e:
        print(f"fid_2={fid_2} zonal统计失败：{e}")
        total_pop = 0

    results.append({"fid_2": fid_2, "pop_sum": total_pop})

# ========== 步骤5：保存Excel ==========
print("正在保存结果为Excel...")
df_result = pd.DataFrame(results)
df_result.to_excel(output_excel_path, index=False)
print(f"已完成：{output_excel_path}")



#计算每个国家不同cdd范围内的人口总数量  次国家尺度 未来2050   人口和cdd 都是1km  arcgis pro  #用的是2030年城镇边界
# -*- coding: utf-8 -*-
import sys
import arcpy
import os
import pandas as pd
import numpy as np

# 设置输出编码为UTF-8，确保可以正确显示中文字符
sys.stdout.reconfigure(encoding='utf-8')

# 设置环境变量
arcpy.env.overwriteOutput = True
arcpy.env.scratchWorkspace = "L:/temp_workspace10"  # 确保路径存在

# 输入文件路径
cdd_raster_path = "M:/2019maxtem/2015-2020data/jisuan18.3_0.1_pop_year_0.125_SUM_change/0ssp126_rasters_1km/re_ssp585_2050_nopop_1km.tif"
pop_raster_path = "L:/worldpop_revise/future_pop_RE_0.125_revise/ssp5/ssp5_RE_2050_0.125_revise_1km.tif" # 更新路径
vector_path = "M:/arcgis_figs3_cityboundary/GHS_UCDB_MTUC_GLOBE_R2024A_GHSL_UCDB_MTUC_2030_GLOBE_R2024_2.shp"

# 输出文件路径
output_excel = "M:/2019maxtem/2015-2020data/jisuan18.3_0.1_pop_year_0.125_SUM_change/0ssp126_rasters_1km/CDD_POP-cities-2050-ssp585-2.xlsx"

# 定义CDD区间，使用一个很大的数代替无穷大
cdd_bins = [0, 600, 1200, 1800, 2400, 3000, 1e10]
cdd_labels = ["lt600", "600_1200", "1200_1800", "1800_2400", "2400_3000", "gte3000"]

# 创建临时工作空间
temp_workspace = "L:/temp_workspace10"
if not os.path.exists(temp_workspace):
    os.makedirs(temp_workspace)
arcpy.env.workspace = temp_workspace

# 重分类CDD栅格到指定区间
reclass_ranges = arcpy.sa.RemapRange([[cdd_bins[i], cdd_bins[i + 1], i + 1] for i in range(len(cdd_bins) - 1)])
try:
    reclass_cdd = arcpy.sa.Reclassify(cdd_raster_path, "VALUE", reclass_ranges, "NODATA")
    reclass_cdd.save(os.path.join(temp_workspace, "reclass_cdd.tif"))
    print("重分类成功")
except arcpy.ExecuteError as e:
    print("重分类时出错: {}".format(e))
    raise

# 以每个国家为单元，对重分类后的CDD栅格和人口栅格进行区域统计
result_df = pd.DataFrame(columns=["Country"] + cdd_labels)

for cdd_label in cdd_labels:
    cdd_zone = cdd_labels.index(cdd_label) + 1
    zone_pop_raster = arcpy.sa.Con(reclass_cdd == cdd_zone, pop_raster_path)
    zonal_table = os.path.join(temp_workspace, "zonal_table_{}.dbf".format(cdd_label))
    try:
        arcpy.sa.ZonalStatisticsAsTable(vector_path, "ID_UC_G0_2", zone_pop_raster, zonal_table, "DATA", "SUM")
    except arcpy.ExecuteError as e:
        print(f"区域统计时出错 (CDD范围: {cdd_label}): {e}")
        continue

    # 读取统计结果并更新DataFrame
    stats_data = arcpy.da.TableToNumPyArray(zonal_table, ["ID_UC_G0_2", "SUM"])
    for row in stats_data:
        country_name = row['ID_UC_G0_2'].strip()
        pop_sum = row['SUM'] if row['SUM'] is not None else 0

        if country_name in result_df["Country"].values:
            result_df.loc[result_df["Country"] == country_name, cdd_label] = pop_sum
        else:
            new_row = pd.DataFrame([{"Country": country_name, cdd_label: pop_sum}])
            result_df = pd.concat([result_df, new_row], ignore_index=True)

# 保存结果到Excel
try:
    result_df.to_excel(output_excel, index=False)
    print("统计完成，结果保存在: {}".format(output_excel))
except Exception as e:
    print(f"保存结果到Excel时出错: {e}")

# 清理临时文件
try:
    arcpy.Delete_management(temp_workspace)
    print("临时文件已清理")
except arcpy.ExecuteError as e:
    print(f"清理临时文件时出错: {e}")




#计算每个国家不同cdd范围内的人口总数量  次国家尺度 未来2050   人口和cdd 都是1km  arcgis pro
#改成用龙赢老师组的1km人口数据
# -*- coding: utf-8 -*-
import sys
import arcpy
import os
import pandas as pd
import numpy as np

# 设置输出编码为UTF-8，确保可以正确显示中文字符
sys.stdout.reconfigure(encoding='utf-8')

# 设置环境变量
arcpy.env.overwriteOutput = True
arcpy.env.scratchWorkspace = "L:/temp_workspace10"  # 确保路径存在

# 输入文件路径
cdd_raster_path = "M:/2019maxtem/2015-2020data/jisuan18.3_0.1_pop_year_0.125_SUM_change/0ssp126_rasters_1km/re_ssp585_2050_nopop_1km.tif"
pop_raster_path = "L:/worldpop_revise/future-1km/SSP5/SSP5_2050.tif" # 更新路径
vector_path = "M:/arcgis_figs3_cityboundary/GHS_UCDB_MTUC_GLOBE_R2024A_GHSL_UCDB_MTUC_2030_GLOBE_R2024_2.shp"

# 输出文件路径
output_excel = "M:/2019maxtem/2015-2020data/jisuan18.3_0.1_pop_year_0.125_SUM_change/0ssp126_rasters_1km/CDD_POP-cities-2050-ssp585-2.xlsx"

# 定义CDD区间，使用一个很大的数代替无穷大
cdd_bins = [0, 600, 1200, 1800, 2400, 3000, 1e10]
cdd_labels = ["lt600", "600_1200", "1200_1800", "1800_2400", "2400_3000", "gte3000"]

# 创建临时工作空间
temp_workspace = "L:/temp_workspace10"
if not os.path.exists(temp_workspace):
    os.makedirs(temp_workspace)
arcpy.env.workspace = temp_workspace

# 重分类CDD栅格到指定区间
reclass_ranges = arcpy.sa.RemapRange([[cdd_bins[i], cdd_bins[i + 1], i + 1] for i in range(len(cdd_bins) - 1)])
try:
    reclass_cdd = arcpy.sa.Reclassify(cdd_raster_path, "VALUE", reclass_ranges, "NODATA")
    reclass_cdd.save(os.path.join(temp_workspace, "reclass_cdd.tif"))
    print("重分类成功")
except arcpy.ExecuteError as e:
    print("重分类时出错: {}".format(e))
    raise

# 以每个国家为单元，对重分类后的CDD栅格和人口栅格进行区域统计
result_df = pd.DataFrame(columns=["Country"] + cdd_labels)

for cdd_label in cdd_labels:
    cdd_zone = cdd_labels.index(cdd_label) + 1
    zone_pop_raster = arcpy.sa.Con(reclass_cdd == cdd_zone, pop_raster_path)
    zonal_table = os.path.join(temp_workspace, "zonal_table_{}.dbf".format(cdd_label))
    try:
        arcpy.sa.ZonalStatisticsAsTable(vector_path, "ID_UC_G0_2", zone_pop_raster, zonal_table, "DATA", "SUM")
    except arcpy.ExecuteError as e:
        print(f"区域统计时出错 (CDD范围: {cdd_label}): {e}")
        continue

    # 读取统计结果并更新DataFrame
    stats_data = arcpy.da.TableToNumPyArray(zonal_table, ["ID_UC_G0_2", "SUM"])
    for row in stats_data:
        country_name = row['ID_UC_G0_2'].strip()
        pop_sum = row['SUM'] if row['SUM'] is not None else 0

        if country_name in result_df["Country"].values:
            result_df.loc[result_df["Country"] == country_name, cdd_label] = pop_sum
        else:
            new_row = pd.DataFrame([{"Country": country_name, cdd_label: pop_sum}])
            result_df = pd.concat([result_df, new_row], ignore_index=True)

# 保存结果到Excel
try:
    result_df.to_excel(output_excel, index=False)
    print("统计完成，结果保存在: {}".format(output_excel))
except Exception as e:
    print(f"保存结果到Excel时出错: {e}")

# 清理临时文件
try:
    arcpy.Delete_management(temp_workspace)
    print("临时文件已清理")
except arcpy.ExecuteError as e:
    print(f"清理临时文件时出错: {e}")






#计算每个国家不同cdd范围内的人口总数量  次国家尺度 未来2050   人口和cdd 都是1km  arcgis pro  #用的是2050年SSP情景，周宇宇团队的
#改成用龙赢老师组的1km人口数据
# -*- coding: utf-8 -*-
import sys
import arcpy
import os
import pandas as pd
import numpy as np

# 设置输出编码为UTF-8，确保可以正确显示中文字符
sys.stdout.reconfigure(encoding='utf-8')

# 设置环境变量
arcpy.env.overwriteOutput = True
arcpy.env.scratchWorkspace = "L:/temp_workspace101"  # 确保路径存在

# 输入文件路径
cdd_raster_path = "M:/2019maxtem/2015-2020data/jisuan18.3_0.1_pop_year_0.125_SUM_change/0ssp126_rasters_1km/re_ssp370_2050_nopop_1km.tif"
pop_raster_path = "L:/worldpop_revise/future-1km/SSP3/SSP3_2050.tif" # 更新路径
vector_path = "M:/ssp城镇化/gUrban_2050_SSP3_SetNull_polygon_编辑.shp"

# 输出文件路径
output_excel = "M:/2019maxtem/2015-2020data/jisuan18.3_0.1_pop_year_0.125_SUM_change/0ssp126_rasters_1km/CDD_POP-cities-2050-ssp370-3.xlsx"

# 定义CDD区间，使用一个很大的数代替无穷大
cdd_bins = [0, 600, 1200, 1800, 2400, 3000, 1e10]
cdd_labels = ["lt600", "600_1200", "1200_1800", "1800_2400", "2400_3000", "gte3000"]

# 创建临时工作空间
temp_workspace = "L:/temp_workspace101"
if not os.path.exists(temp_workspace):
    os.makedirs(temp_workspace)
arcpy.env.workspace = temp_workspace

# 重分类CDD栅格到指定区间
reclass_ranges = arcpy.sa.RemapRange([[cdd_bins[i], cdd_bins[i + 1], i + 1] for i in range(len(cdd_bins) - 1)])
try:
    reclass_cdd = arcpy.sa.Reclassify(cdd_raster_path, "VALUE", reclass_ranges, "NODATA")
    reclass_cdd.save(os.path.join(temp_workspace, "reclass_cdd.tif"))
    print("重分类成功")
except arcpy.ExecuteError as e:
    print("重分类时出错: {}".format(e))
    raise

# 以每个国家为单元，对重分类后的CDD栅格和人口栅格进行区域统计
result_df = pd.DataFrame(columns=["Country"] + cdd_labels)

for cdd_label in cdd_labels:
    cdd_zone = cdd_labels.index(cdd_label) + 1
    zone_pop_raster = arcpy.sa.Con(reclass_cdd == cdd_zone, pop_raster_path)
    zonal_table = os.path.join(temp_workspace, "zonal_table_{}.dbf".format(cdd_label))
    try:
        arcpy.sa.ZonalStatisticsAsTable(vector_path, "name", zone_pop_raster, zonal_table, "DATA", "SUM")
    except arcpy.ExecuteError as e:
        print(f"区域统计时出错 (CDD范围: {cdd_label}): {e}")
        continue

    # 读取统计结果并更新DataFrame
    stats_data = arcpy.da.TableToNumPyArray(zonal_table, ["name", "SUM"])
    for row in stats_data:
        country_name = row['name'].strip()
        pop_sum = row['SUM'] if row['SUM'] is not None else 0

        if country_name in result_df["Country"].values:
            result_df.loc[result_df["Country"] == country_name, cdd_label] = pop_sum
        else:
            new_row = pd.DataFrame([{"Country": country_name, cdd_label: pop_sum}])
            result_df = pd.concat([result_df, new_row], ignore_index=True)

# 保存结果到Excel
try:
    result_df.to_excel(output_excel, index=False)
    print("统计完成，结果保存在: {}".format(output_excel))
except Exception as e:
    print(f"保存结果到Excel时出错: {e}")

# 清理临时文件
try:
    arcpy.Delete_management(temp_workspace)
    print("临时文件已清理")
except arcpy.ExecuteError as e:
    print(f"清理临时文件时出错: {e}")





#计算每个国家不同cdd范围内的人口总数量  次国家尺度 未来2050   人口和cdd 都是1km  arcgis pro  #用的是2050年SSP情景，黎夏团队的
#改成用龙赢老师组的1km人口数据
# -*- coding: utf-8 -*-
import sys
import arcpy
import os
import pandas as pd
import numpy as np

# 设置输出编码为UTF-8，确保可以正确显示中文字符
sys.stdout.reconfigure(encoding='utf-8')

# 设置环境变量
arcpy.env.overwriteOutput = True
arcpy.env.scratchWorkspace = "L:/temp_workspace102"  # 确保路径存在

# 输入文件路径
cdd_raster_path = "M:/2019maxtem/2015-2020data/jisuan18.3_0.1_pop_year_0.125_SUM_change/0ssp126_rasters_1km/re_ssp585_2050_nopop_1km.tif"
pop_raster_path = "L:/worldpop_revise/future-1km/SSP5/SSP5_2050.tif" # 更新路径
vector_path = "M:/ssp城镇化-lixia/gUrban_2050_SSP5_SetNull_polygon1_over1km_identity.shp"

# 输出文件路径
output_excel = "M:/2019maxtem/2015-2020data/jisuan18.3_0.1_pop_year_0.125_SUM_change/0ssp126_rasters_1km/CDD_POP-cities-2050-ssp585-4.xlsx"

# 定义CDD区间，使用一个很大的数代替无穷大
cdd_bins = [0, 600, 1200, 1800, 2400, 3000, 1e10]
cdd_labels = ["lt600", "600_1200", "1200_1800", "1800_2400", "2400_3000", "gte3000"]

# 创建临时工作空间
temp_workspace = "L:/temp_workspace102"
if not os.path.exists(temp_workspace):
    os.makedirs(temp_workspace)
arcpy.env.workspace = temp_workspace

# 重分类CDD栅格到指定区间
reclass_ranges = arcpy.sa.RemapRange([[cdd_bins[i], cdd_bins[i + 1], i + 1] for i in range(len(cdd_bins) - 1)])
try:
    reclass_cdd = arcpy.sa.Reclassify(cdd_raster_path, "VALUE", reclass_ranges, "NODATA")
    reclass_cdd.save(os.path.join(temp_workspace, "reclass_cdd.tif"))
    print("重分类成功")
except arcpy.ExecuteError as e:
    print("重分类时出错: {}".format(e))
    raise

# 以每个国家为单元，对重分类后的CDD栅格和人口栅格进行区域统计
result_df = pd.DataFrame(columns=["Country"] + cdd_labels)

for cdd_label in cdd_labels:
    cdd_zone = cdd_labels.index(cdd_label) + 1
    zone_pop_raster = arcpy.sa.Con(reclass_cdd == cdd_zone, pop_raster_path)
    zonal_table = os.path.join(temp_workspace, "zonal_table_{}.dbf".format(cdd_label))
    try:
        arcpy.sa.ZonalStatisticsAsTable(vector_path, "name", zone_pop_raster, zonal_table, "DATA", "SUM")
    except arcpy.ExecuteError as e:
        print(f"区域统计时出错 (CDD范围: {cdd_label}): {e}")
        continue

    # 读取统计结果并更新DataFrame
    stats_data = arcpy.da.TableToNumPyArray(zonal_table, ["name", "SUM"])
    for row in stats_data:
        country_name = row['name'].strip()
        pop_sum = row['SUM'] if row['SUM'] is not None else 0

        if country_name in result_df["Country"].values:
            result_df.loc[result_df["Country"] == country_name, cdd_label] = pop_sum
        else:
            new_row = pd.DataFrame([{"Country": country_name, cdd_label: pop_sum}])
            result_df = pd.concat([result_df, new_row], ignore_index=True)

# 保存结果到Excel
try:
    result_df.to_excel(output_excel, index=False)
    print("统计完成，结果保存在: {}".format(output_excel))
except Exception as e:
    print(f"保存结果到Excel时出错: {e}")

# 清理临时文件
try:
    arcpy.Delete_management(temp_workspace)
    print("临时文件已清理")
except arcpy.ExecuteError as e:
    print(f"清理临时文件时出错: {e}")





#尝试
# -*- coding: utf-8 -*-
import sys
import arcpy
import os
import pandas as pd
import numpy as np

# 设置输出编码为UTF-8，确保可以正确显示中文字符
sys.stdout.reconfigure(encoding='utf-8')

# 设置环境变量
arcpy.env.overwriteOutput = True
arcpy.env.scratchWorkspace = "L:/temp_workspace10"  # 确保路径存在

# 输入文件路径
cdd_raster_path = "M:/2019maxtem/2015-2020data/jisuan18.3_0.1_pop_year_0.125_SUM_change/0ssp126_rasters_1km/re_ssp370_2050_nopop_1km.tif"
pop_raster_path = "L:/worldpop_revise/future_pop_RE_0.125_revise/ssp3/ssp3_RE_2050_0.125_revise_1km.tif" # 更新路径
vector_path = "M:/贫困-最新发布/GSAP_AM24_2019/GSAP_AM24_2019.shp"

# 输出文件路径
output_excel = "M:/2019maxtem/2015-2020data/jisuan18.3_0.1_pop_year_0.125_SUM_change/0ssp126_rasters_1km/CDD_POP-ciguojia-2050-ssp370-changshi.xlsx"

# 定义CDD区间，使用一个很大的数代替无穷大
cdd_bins = [0, 600, 1200, 1800, 2400, 3000, 1e10]
cdd_labels = ["lt600", "600_1200", "1200_1800", "1800_2400", "2400_3000", "gte3000"]

# 创建临时工作空间
temp_workspace = "L:/temp_workspace10"
if not os.path.exists(temp_workspace):
    os.makedirs(temp_workspace)
arcpy.env.workspace = temp_workspace

# 重分类CDD栅格到指定区间
reclass_ranges = arcpy.sa.RemapRange([[cdd_bins[i], cdd_bins[i + 1], i + 1] for i in range(len(cdd_bins) - 1)])
try:
    reclass_cdd = arcpy.sa.Reclassify(cdd_raster_path, "VALUE", reclass_ranges, "NODATA")
    reclass_cdd.save(os.path.join(temp_workspace, "reclass_cdd.tif"))
    print("重分类成功")
except arcpy.ExecuteError as e:
    print("重分类时出错: {}".format(e))
    raise

# 以每个国家为单元，对重分类后的CDD栅格和人口栅格进行区域统计
result_df = pd.DataFrame(columns=["Country"] + cdd_labels)

for cdd_label in cdd_labels:
    cdd_zone = cdd_labels.index(cdd_label) + 1
    zone_pop_raster = arcpy.sa.Con(reclass_cdd == cdd_zone, pop_raster_path)
    zonal_table = os.path.join(temp_workspace, "zonal_table_{}.dbf".format(cdd_label))
    try:
        arcpy.sa.ZonalStatisticsAsTable(vector_path, "geo_code", zone_pop_raster, zonal_table, "DATA", "SUM")
    except arcpy.ExecuteError as e:
        print(f"区域统计时出错 (CDD范围: {cdd_label}): {e}")
        continue

    # 读取统计结果并更新DataFrame
    stats_data = arcpy.da.TableToNumPyArray(zonal_table, ["geo_code", "SUM"])
    for row in stats_data:
        country_name = row['geo_code'].strip()
        pop_sum = row['SUM'] if row['SUM'] is not None else 0

        if country_name in result_df["Country"].values:
            result_df.loc[result_df["Country"] == country_name, cdd_label] = pop_sum
        else:
            new_row = pd.DataFrame([{"Country": country_name, cdd_label: pop_sum}])
            result_df = pd.concat([result_df, new_row], ignore_index=True)

# 保存结果到Excel
try:
    result_df.to_excel(output_excel, index=False)
    print("统计完成，结果保存在: {}".format(output_excel))
except Exception as e:
    print(f"保存结果到Excel时出错: {e}")

# 清理临时文件
try:
    arcpy.Delete_management(temp_workspace)
    print("临时文件已清理")
except arcpy.ExecuteError as e:
    print(f"清理临时文件时出错: {e}")






#对未来2050年人口数据进行降尺度处理  根据2020年数据进行降尺度，是可用的  不过缺少对结果的检验
# -*- coding: utf-8 -*-
import arcpy
from arcpy.sa import *
import os
# 初始化 ArcGIS 环境
arcpy.env.overwriteOutput = True
arcpy.CheckOutExtension("Spatial")
# === 输入路径 ===
hist_1km_path = r"F:\人口数据\lanscan_origin\landscan-global-2021-assets\landscan-global-2021_360.tif"
future_0125_path = r"L:\worldpop_revise\future_pop_RE_0.125_revise\ssp5\ssp5_RE_2050_0.125_revise.tif"
# === 输出路径 ===
output_downscaled = r"L:\worldpop_revise\future_pop_RE_0.125_revise\ssp5\ssp5_RE_2050_0.125_revise_1km.tif"
temp_workspace = r"F:\temp_downscale"
if not os.path.exists(temp_workspace):
    os.makedirs(temp_workspace)
arcpy.env.workspace = temp_workspace
# === 步骤 1：读取栅格 ===
hist_raster = Raster(hist_1km_path)
fut_raster = Raster(future_0125_path)
# === 步骤 2：汇总历史人口至0.125°（15x15 像元） ===
agg_factor = 15
agg_hist = BlockStatistics(hist_raster, NbrRectangle(agg_factor, agg_factor, "CELL"), "SUM")
agg_hist.save("agg_hist")
# === 步骤 3：计算比例并去除无穷大或异常值 ===
raw_ratio = Divide(hist_raster, agg_hist + 1e-6)
hist_ratio = Con(raw_ratio > 100, 0, raw_ratio)  # 超过100的人口比例设为0
hist_ratio.save(r"M:\hist_ratio_clean.tif")
# === 步骤 4：上采样未来人口数据到 1km ===
fut_resampled = arcpy.management.Resample(fut_raster, "fut_resample", hist_raster.meanCellWidth, "NEAREST")[0]
# === 步骤 5：按比例分配未来人口 ===
downscaled_raster = Times(fut_resampled, hist_ratio)
downscaled_raster.save(output_downscaled)
print("✅ 降尺度完成，结果保存至：", output_downscaled)





#对未来2050年cdd数据进行降尺度处理  根据2020年数据进行降尺度，是可用的  不过缺少对结果的检验
# -*- coding: utf-8 -*-
import arcpy
from arcpy.sa import *
import os

# 初始化 ArcGIS 环境
arcpy.env.overwriteOutput = True
arcpy.CheckOutExtension("Spatial")

# === 输入路径 ===
hist_1km_path = r"M:\2019maxtem\2015-2020data\jisuan18.3_0.1_pop_year_0.125_SUM_change\0ssp126_rasters_1km\re_ssp126_2021.tif"
future_0125_path = r"M:\2019maxtem\2015-2020data\jisuan18.3_0.1_pop_year_0.125_SUM_change\0ssp126_rasters_1km\re_ssp585_2050_nopop.tif"

# === 输出路径 ===
output_downscaled = r"M:\2019maxtem\2015-2020data\jisuan18.3_0.1_pop_year_0.125_SUM_change\0ssp126_rasters_1km\re_ssp585_2050_nopop_1km.tif"
temp_workspace = r"M:\temp_downscale_cdd"
if not os.path.exists(temp_workspace):
    os.makedirs(temp_workspace)
arcpy.env.workspace = temp_workspace

# === 步骤 1：读取历史与未来栅格 ===
hist_raster = Raster(hist_1km_path)
fut_raster = Raster(future_0125_path)

# === 步骤 2：将未来数据从平均值转换为总量（乘以 225）===
fut_total = Times(fut_raster, 225)

# === 步骤 3：聚合历史数据（15×15）作为大格总量参考 ===
agg_hist = BlockStatistics(hist_raster, NbrRectangle(15, 15, "CELL"), "SUM")
agg_hist.save(os.path.join(temp_workspace, "agg_hist.tif"))

# === 步骤 4：将聚合结果重采样回 1km，构建比例因子 ===
agg_hist_resample = arcpy.management.Resample(
    agg_hist,
    os.path.join(temp_workspace, "agg_hist_resample.tif"),
    hist_raster.meanCellWidth,
    "NEAREST"
)[0]
agg_hist_resample = Raster(agg_hist_resample)

# === 步骤 5：构建比例图，处理异常值、避免除以零 ===
ratio_raw = Divide(hist_raster, agg_hist_resample + 1e-6)
hist_ratio = Con((ratio_raw > 0) & (ratio_raw < 100), ratio_raw, 0)
hist_ratio.save(os.path.join(temp_workspace, "hist_ratio_clean.tif"))

# === 步骤 6：将未来总量上采样为 1km 分辨率 ===
fut_resampled = arcpy.management.Resample(
    fut_total,
    os.path.join(temp_workspace, "fut_total_resample.tif"),
    hist_raster.meanCellWidth,
    "NEAREST"
)[0]
fut_resampled = Raster(fut_resampled)

# === 步骤 7：按比例进行降尺度分配（总量 × 比例） ===
final_downscaled = Times(fut_resampled, hist_ratio)
final_downscaled.save(output_downscaled)

print("\n✅ 降尺度完成，结果保存至：", output_downscaled)

