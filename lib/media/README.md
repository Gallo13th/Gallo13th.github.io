# `public/lib/media/` —— 视听页播放器的自托管第三方库

这三个文件是**手工落库**的第三方产物（不是构建时从 `node_modules/` 生成的），
因此放在 `public/lib/` 而不是 `public/vendor/` —— 后者是 `.gitignore` 的、由 `prebuild` 生成。

## 为什么不用 CDN

原先 `/media/` 的播放器从 `cdn.jsdelivr.net` 载入。实测 jsDelivr 对 `@main` 之类的引用会
**301 到 `raw.githubusercontent.com`**，而该域在中国大陆常见不可达 —— 即读者/作者很可能点了
「加载播放器」却什么也拿不到。这三个文件是页面唯一还会走外部 CDN 的资源，改为自托管后，
`/media/` 的首屏与播放器都不再依赖任何 CDN。

## 来源与校验

| 文件 | 版本 | 上游（固定版本 URL） | 字节 | SHA-256 |
|---|---|---|---|---|
| `APlayer.min.css` | aplayer@1.10.1 | `https://cdn.jsdelivr.net/npm/aplayer@1.10.1/dist/APlayer.min.css` | 12528 | `baa4101a70dc9912af84ac1ce559b85d3d46436a15eadd54d0d47637db55f814` |
| `APlayer.min.js` | aplayer@1.10.1 | `https://cdn.jsdelivr.net/npm/aplayer@1.10.1/dist/APlayer.min.js` | 59325 | `e98ec22436a5b6878d824f997ed8020fd8cb8261afe31294a3c9d0d07800c15a` |
| `Meting.min.js` | meting@2.0.2 | `https://cdn.jsdelivr.net/npm/meting@2.0.2/dist/Meting.min.js` | 2908 | `eb9e8f9f495fcdb8c583313a39db04905b1ac65d327b81d93b357c7c7f3f9d70` |

许可证：**两者均为 MIT**（`aplayer` 与 `meting` 的 npm `package.json` 中 `license: MIT`）。
上游仓库：<https://github.com/DIYgod/APlayer>、<https://github.com/metowolf/MetingJS>。
（npm registry 在本机不可达，故未走 `npm install` 路线；文件已入库并记下哈希，
`scripts/check-media-player.mjs` 会校验它们存在且非空。）

## 两个已核验的性质（别再猜）

1. **`APlayer.min.css` 不含任何 `url(...)`** —— 图标是 JS 内联 SVG，没有需要一并搬运的字体/图片。
   `check-media-player.mjs` 把这条锁成了断言：将来升级 APlayer 若开始引用外部资源，门会红。
2. **音频字节最终来自网易云自己的 CDN**：数据 API 返回的 `url` 是一个 **302**，
   `location` 指向 `m7.music.126.net`。也就是说第三方实例只负责
   **歌单元信息 + 带签名的地址**（实测 871 B / 2 首歌），音频流不经过它。

## 仍然存在的第三方依赖（**删不掉，已改为可见**）

播放器的**数据**来自 Meting 的公共实例 `https://api.i-meto.com/meting/api`。
自托管 JS 不能消除它 —— 除非自建一个 Meting API 服务端，或改用网易云官方外链播放器
（实测 `outchain/player?type=0|2|3` 对匿名请求返回**完全相同的 5673 B 通用外壳**，
无法确认歌单能否正常播放，故未采用）。

因此 `/media/` 的加载流程里加了一次**预检**：先 `fetch` 这个 API（实测带
`access-control-allow-origin: *`，可跨域预检），拿不到或返回空列表时**明确提示并给出
「在网易云打开」的出路**，而不是让 Meting 安静地渲染一个空播放器。
`media/index.astro` 里那个 URL 只出现一次（`METING_API` 常量），要换服务端只改一处。
